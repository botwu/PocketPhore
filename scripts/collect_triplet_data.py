"""
蛋白-配体-药效团三元组数据收集脚本

数据流程:
1. 下载PDB复合物
2. 提取蛋白口袋残基
3. 生成药效团模型
4. 保存为三元组格式
"""

import os
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np
from tqdm import tqdm
import pickle
import json

class PocketSelector(Select):
    """选择配体周围的口袋残基"""
    def __init__(self, ligand_coords, cutoff=10.0):
        self.ligand_coords = ligand_coords
        self.cutoff = cutoff

    def accept_residue(self, residue):
        # 计算残基到配体的最小距离
        min_dist = float('inf')
        for atom in residue:
            atom_coord = atom.get_coord()
            for lig_coord in self.ligand_coords:
                dist = np.linalg.norm(atom_coord - lig_coord)
                min_dist = min(min_dist, dist)

        return min_dist < self.cutoff

def download_pdb_complex(pdb_id, output_dir='./data/pdb_raw'):
    """下载PDB复合物文件"""
    os.makedirs(output_dir, exist_ok=True)

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_path = os.path.join(output_dir, f"{pdb_id}.pdb")

    if os.path.exists(output_path):
        return output_path

    try:
        response = requests.get(url)
        response.raise_for_status()

        with open(output_path, 'w') as f:
            f.write(response.text)

        return output_path
    except Exception as e:
        print(f"下载 {pdb_id} 失败: {e}")
        return None

def extract_ligand_from_pdb(pdb_file, ligand_resname='LIG'):
    """从PDB文件提取配体坐标和SMILES"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    ligand_coords = []
    ligand_atoms = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == ligand_resname:
                    for atom in residue:
                        ligand_coords.append(atom.get_coord())
                        ligand_atoms.append(atom.element)

    return np.array(ligand_coords), ligand_atoms

def extract_pocket_residues(pdb_file, ligand_coords, cutoff=10.0,
                           output_dir='./data/pockets'):
    """提取配体周围的口袋残基"""
    os.makedirs(output_dir, exist_ok=True)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    # 选择口袋残基
    io = PDBIO()
    io.set_structure(structure)

    pocket_selector = PocketSelector(ligand_coords, cutoff)

    pdb_id = os.path.basename(pdb_file).replace('.pdb', '')
    pocket_file = os.path.join(output_dir, f"{pdb_id}_pocket.pdb")

    io.save(pocket_file, pocket_selector)

    return pocket_file

def generate_pharmacophore_from_ligand(ligand_coords, ligand_mol,
                                      output_dir='./data/pharmacophores'):
    """
    从配体生成药效团模型

    方法1: 使用RDKit自动检测药效团特征
    方法2: 使用LigandScout (如果有license)
    方法3: 使用开源工具Pharao
    """
    os.makedirs(output_dir, exist_ok=True)

    from rdkit.Chem import ChemicalFeatures
    from rdkit import RDConfig

    # 加载药效团特征定义
    fdef_name = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdef_name)

    # 检测药效团特征
    features = factory.GetFeaturesForMol(ligand_mol)

    pharmacophore = []
    for feat in features:
        feat_type = feat.GetFamily()
        feat_pos = feat.GetPos()

        # 映射RDKit特征到PhoreGen特征类型
        phore_type_map = {
            'Donor': 'HD',           # 氢键供体
            'Acceptor': 'HA',        # 氢键受体
            'Aromatic': 'AR',        # 芳香环
            'Hydrophobe': 'HY',      # 疏水
            'PosIonizable': 'PO',    # 正电
            'NegIonizable': 'NE'     # 负电
        }

        if feat_type in phore_type_map:
            pharmacophore.append({
                'type': phore_type_map[feat_type],
                'position': [feat_pos.x, feat_pos.y, feat_pos.z],
                'radius': 1.0  # 默认容忍半径
            })

    return pharmacophore

def process_pdbbind_entry(pdb_id, pdbbind_dir='./data/pdbbind'):
    """
    处理PDBbind数据集中的一个条目

    PDBbind目录结构:
    pdbbind/
        v2020/
            {pdb_id}/
                {pdb_id}_protein.pdb
                {pdb_id}_ligand.mol2
                {pdb_id}_ligand.sdf
    """
    try:
        # 1. 读取配体
        ligand_file = os.path.join(pdbbind_dir, 'v2020', pdb_id,
                                   f'{pdb_id}_ligand.sdf')

        if not os.path.exists(ligand_file):
            print(f"配体文件不存在: {ligand_file}")
            return None

        ligand_mol = Chem.SDMolSupplier(ligand_file)[0]
        if ligand_mol is None:
            return None

        # 获取配体3D坐标
        conf = ligand_mol.GetConformer()
        ligand_coords = conf.GetPositions()

        # 2. 读取蛋白质
        protein_file = os.path.join(pdbbind_dir, 'v2020', pdb_id,
                                   f'{pdb_id}_protein.pdb')

        if not os.path.exists(protein_file):
            return None

        # 3. 提取口袋残基
        pocket_file = extract_pocket_residues(protein_file, ligand_coords,
                                             cutoff=10.0)

        # 4. 生成药效团
        pharmacophore = generate_pharmacophore_from_ligand(ligand_coords,
                                                          ligand_mol)

        # 5. 构建三元组
        triplet = {
            'pdb_id': pdb_id,
            'ligand': {
                'mol': ligand_mol,
                'smiles': Chem.MolToSmiles(ligand_mol),
                'coords': ligand_coords.tolist(),
                'file': ligand_file
            },
            'protein': {
                'pocket_file': pocket_file,
                'full_file': protein_file
            },
            'pharmacophore': pharmacophore
        }

        return triplet

    except Exception as e:
        print(f"处理 {pdb_id} 时出错: {e}")
        return None

def collect_pdbbind_triplets(pdbbind_dir='./data/pdbbind',
                            output_file='./data/triplets.pkl'):
    """
    批量处理PDBbind数据集
    """
    # 读取PDBbind索引文件
    index_file = os.path.join(pdbbind_dir, 'v2020', 'index',
                             'INDEX_refined_data.2020')

    pdb_ids = []
    with open(index_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) >= 1:
                pdb_ids.append(parts[0])

    print(f"找到 {len(pdb_ids)} 个PDB条目")

    # 处理每个条目
    triplets = []
    for pdb_id in tqdm(pdb_ids, desc="处理PDBbind数据"):
        triplet = process_pdbbind_entry(pdb_id, pdbbind_dir)
        if triplet is not None:
            triplets.append(triplet)

    print(f"成功收集 {len(triplets)} 个三元组")

    # 保存结果
    with open(output_file, 'wb') as f:
        pickle.dump(triplets, f)

    print(f"数据已保存至: {output_file}")

    return triplets

def split_dataset(triplets, train_ratio=0.8, val_ratio=0.1):
    """划分训练/验证/测试集"""
    n = len(triplets)
    indices = np.random.permutation(n)

    n_train = int(n * train_ratio)
    n_val = int(n * val_ratio)

    train_indices = indices[:n_train]
    val_indices = indices[n_train:n_train+n_val]
    test_indices = indices[n_train+n_val:]

    splits = {
        'train': [triplets[i] for i in train_indices],
        'val': [triplets[i] for i in val_indices],
        'test': [triplets[i] for i in test_indices]
    }

    return splits

if __name__ == "__main__":
    # 方案1: 使用PDBbind数据集 (推荐)
    print("开始收集蛋白-配体-药效团三元组数据...")

    # 需要先下载PDBbind数据集
    # http://www.pdbbind.org.cn/download.php
    pdbbind_dir = './data/pdbbind'

    if os.path.exists(pdbbind_dir):
        triplets = collect_pdbbind_triplets(pdbbind_dir)

        # 划分数据集
        splits = split_dataset(triplets)

        print(f"训练集: {len(splits['train'])} 条")
        print(f"验证集: {len(splits['val'])} 条")
        print(f"测试集: {len(splits['test'])} 条")
    else:
        print(f"请先下载PDBbind数据集到 {pdbbind_dir}")
        print("下载地址: http://www.pdbbind.org.cn/download.php")
