"""
完整的蛋白-配体-药效团三元组数据收集Pipeline

使用方法:
python full_pipeline.py --pdbbind_dir ./data/pdbbind --output_dir ./data/processed

流程:
1. 从PDBbind加载蛋白-配体复合物
2. 提取口袋残基
3. 生成药效团
4. 提取口袋特征
5. 质量过滤
6. 保存最终数据集
"""

import os
import argparse
import pickle
import numpy as np
from tqdm import tqdm
import multiprocessing as mp
from rdkit import Chem

# 导入前面的模块
import sys
sys.path.append('./scripts')
from extract_pocket_features import PocketFeatureExtractor
from filter_triplet_data import TripletQualityFilter

def process_single_entry(args):
    """处理单个PDBbind条目 (用于并行处理)"""
    pdb_id, pdbbind_dir, use_ligandscout = args

    try:
        # 1. 读取配体
        ligand_file = os.path.join(pdbbind_dir, 'v2020', pdb_id,
                                   f'{pdb_id}_ligand.sdf')

        if not os.path.exists(ligand_file):
            return None

        ligand_suppl = Chem.SDMolSupplier(ligand_file)
        ligand_mol = ligand_suppl[0]

        if ligand_mol is None:
            return None

        # 获取配体坐标
        conf = ligand_mol.GetConformer()
        ligand_coords = conf.GetPositions()

        # 2. 读取蛋白质
        protein_file = os.path.join(pdbbind_dir, 'v2020', pdb_id,
                                   f'{pdb_id}_protein.pdb')

        if not os.path.exists(protein_file):
            return None

        # 3. 生成药效团
        if use_ligandscout:
            # 使用LigandScout (如果安装)
            from pharmacophore_ligandscout import generate_pharmacophore_ligandscout
            pharmacophore = generate_pharmacophore_ligandscout(
                protein_file, ligand_file,
                f'./data/pharmacophores/{pdb_id}.pml'
            )
        else:
            # 使用RDKit (默认)
            pharmacophore = generate_pharmacophore_rdkit(ligand_mol)

        if not pharmacophore or len(pharmacophore) == 0:
            return None

        # 4. 提取口袋特征
        extractor = PocketFeatureExtractor(protein_file, ligand_coords,
                                          pocket_cutoff=10.0)
        pocket_features = extractor.extract_all_features()

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
                'pocket_file': protein_file,
            },
            'pharmacophore': pharmacophore,
            'pocket_features': pocket_features
        }

        return triplet

    except Exception as e:
        print(f"处理 {pdb_id} 时出错: {e}")
        return None

def generate_pharmacophore_rdkit(ligand_mol):
    """使用RDKit生成药效团 (简化版)"""
    from rdkit.Chem import ChemicalFeatures
    from rdkit import RDConfig

    fdef_name = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdef_name)

    features = factory.GetFeaturesForMol(ligand_mol)

    pharmacophore = []

    phore_type_map = {
        'Donor': 'HD',
        'Acceptor': 'HA',
        'Aromatic': 'AR',
        'Hydrophobe': 'HY',
        'PosIonizable': 'PO',
        'NegIonizable': 'NE'
    }

    for feat in features:
        feat_type = feat.GetFamily()
        feat_pos = feat.GetPos()

        if feat_type in phore_type_map:
            pharmacophore.append({
                'type': phore_type_map[feat_type],
                'position': [feat_pos.x, feat_pos.y, feat_pos.z],
                'radius': 1.0
            })

    return pharmacophore

def collect_pdbbind_parallel(pdbbind_dir, output_dir, n_workers=8,
                            use_ligandscout=False):
    """并行处理PDBbind数据集"""
    os.makedirs(output_dir, exist_ok=True)

    # 读取PDBbind索引
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

    # 并行处理
    args_list = [(pdb_id, pdbbind_dir, use_ligandscout) for pdb_id in pdb_ids]

    with mp.Pool(n_workers) as pool:
        results = list(tqdm(
            pool.imap(process_single_entry, args_list),
            total=len(args_list),
            desc="处理PDBbind数据"
        ))

    # 过滤None结果
    triplets = [r for r in results if r is not None]

    print(f"成功处理 {len(triplets)} 个三元组")

    # 保存原始数据
    raw_file = os.path.join(output_dir, 'triplets_raw.pkl')
    with open(raw_file, 'wb') as f:
        pickle.dump(triplets, f)

    print(f"原始数据已保存至: {raw_file}")

    return triplets

def main():
    parser = argparse.ArgumentParser(description='完整数据收集Pipeline')
    parser.add_argument('--pdbbind_dir', type=str, required=True,
                       help='PDBbind数据集路径')
    parser.add_argument('--output_dir', type=str, default='./data/processed',
                       help='输出目录')
    parser.add_argument('--n_workers', type=int, default=8,
                       help='并行进程数')
    parser.add_argument('--use_ligandscout', action='store_true',
                       help='使用LigandScout生成药效团')
    parser.add_argument('--skip_collection', action='store_true',
                       help='跳过数据收集，直接进行过滤')

    args = parser.parse_args()

    print("="*60)
    print("蛋白-配体-药效团三元组数据收集Pipeline")
    print("="*60)

    # Step 1: 数据收集
    if not args.skip_collection:
        print("\n[1/3] 数据收集...")
        triplets = collect_pdbbind_parallel(
            args.pdbbind_dir,
            args.output_dir,
            n_workers=args.n_workers,
            use_ligandscout=args.use_ligandscout
        )
    else:
        print("\n[1/3] 跳过数据收集，加载现有数据...")
        raw_file = os.path.join(args.output_dir, 'triplets_raw.pkl')
        with open(raw_file, 'rb') as f:
            triplets = pickle.load(f)

    # Step 2: 质量过滤
    print("\n[2/3] 质量过滤...")
    filter = TripletQualityFilter()
    filtered_triplets, stats = filter.filter_dataset(triplets)

    # 保存过滤后的数据
    filtered_file = os.path.join(args.output_dir, 'triplets_filtered.pkl')
    with open(filtered_file, 'wb') as f:
        pickle.dump(filtered_triplets, f)

    print(f"过滤后数据已保存至: {filtered_file}")

    # Step 3: 数据集划分
    print("\n[3/3] 数据集划分...")
    n = len(filtered_triplets)
    indices = np.random.permutation(n)

    n_train = int(n * 0.8)
    n_val = int(n * 0.1)

    splits = {
        'train': [filtered_triplets[i] for i in indices[:n_train]],
        'val': [filtered_triplets[i] for i in indices[n_train:n_train+n_val]],
        'test': [filtered_triplets[i] for i in indices[n_train+n_val:]]
    }

    # 保存划分后的数据
    for split_name, split_data in splits.items():
        split_file = os.path.join(args.output_dir, f'{split_name}.pkl')
        with open(split_file, 'wb') as f:
            pickle.dump(split_data, f)

        print(f"{split_name}集: {len(split_data)} 条 -> {split_file}")

    print("\n" + "="*60)
    print("数据收集完成!")
    print("="*60)
    print(f"总数据量: {n}")
    print(f"训练集: {len(splits['train'])}")
    print(f"验证集: {len(splits['val'])}")
    print(f"测试集: {len(splits['test'])}")
    print("="*60)

if __name__ == "__main__":
    main()
