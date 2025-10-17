"""
蛋白质口袋残基特征提取模块

提取特征:
1. 残基类型
2. 二级结构 (使用DSSP)
3. 疏水性/极性
4. 侧链朝向
5. 几何特征 (法向量、曲率)
"""

import numpy as np
from Bio.PDB import PDBParser, DSSP, Selection
from Bio.PDB.vectors import Vector, calc_angle
from scipy.spatial import KDTree
import warnings
warnings.filterwarnings('ignore')

# 氨基酸类型到索引的映射
AA_TO_INDEX = {
    'ALA': 0, 'CYS': 1, 'ASP': 2, 'GLU': 3, 'PHE': 4,
    'GLY': 5, 'HIS': 6, 'ILE': 7, 'LYS': 8, 'LEU': 9,
    'MET': 10, 'ASN': 11, 'PRO': 12, 'GLN': 13, 'ARG': 14,
    'SER': 15, 'THR': 16, 'VAL': 17, 'TRP': 18, 'TYR': 19,
    'UNK': 20  # 未知
}

# Kyte-Doolittle疏水性量表
HYDROPHOBICITY_SCALE = {
    'ALA': 1.8, 'CYS': 2.5, 'ASP': -3.5, 'GLU': -3.5, 'PHE': 2.8,
    'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5, 'LYS': -3.9, 'LEU': 3.8,
    'MET': 1.9, 'ASN': -3.5, 'PRO': -1.6, 'GLN': -3.5, 'ARG': -4.5,
    'SER': -0.8, 'THR': -0.7, 'VAL': 4.2, 'TRP': -0.9, 'TYR': -1.3
}

# 二级结构映射
SS_TO_INDEX = {
    'H': 0,  # α-螺旋
    'B': 1,  # β-折叠
    'E': 1,  # β-折叠延伸
    'G': 0,  # 3-10螺旋
    'I': 0,  # π-螺旋
    'T': 2,  # Turn
    'S': 2,  # Bend
    '-': 2,  # 无规卷曲
    ' ': 2   # 未定义
}

class PocketFeatureExtractor:
    def __init__(self, pdb_file, ligand_coords, pocket_cutoff=10.0):
        """
        初始化口袋特征提取器

        Args:
            pdb_file: PDB文件路径
            ligand_coords: 配体坐标 (N, 3)
            pocket_cutoff: 口袋半径 (Å)
        """
        self.pdb_file = pdb_file
        self.ligand_coords = ligand_coords
        self.pocket_cutoff = pocket_cutoff

        # 解析PDB结构
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure('protein', pdb_file)

        # 提取口袋残基
        self.pocket_residues = self._extract_pocket_residues()

        # 计算二级结构
        self.dssp_data = self._compute_dssp()

        # 计算口袋中心
        self.pocket_center = self._compute_pocket_center()

    def _extract_pocket_residues(self):
        """提取配体周围的口袋残基"""
        pocket_residues = []

        for model in self.structure:
            for chain in model:
                for residue in chain:
                    # 跳过水分子和配体
                    if residue.id[0] != ' ':
                        continue

                    # 计算残基到配体的最小距离
                    min_dist = float('inf')
                    for atom in residue:
                        atom_coord = atom.get_coord()
                        for lig_coord in self.ligand_coords:
                            dist = np.linalg.norm(atom_coord - lig_coord)
                            min_dist = min(min_dist, dist)

                    if min_dist < self.pocket_cutoff:
                        pocket_residues.append(residue)

        print(f"找到 {len(pocket_residues)} 个口袋残基")
        return pocket_residues

    def _compute_dssp(self):
        """使用DSSP计算二级结构"""
        try:
            model = self.structure[0]
            dssp = DSSP(model, self.pdb_file, dssp='mkdssp')
            return dssp
        except Exception as e:
            print(f"DSSP计算失败: {e}")
            print("请安装DSSP: conda install -c salilab dssp")
            return None

    def _compute_pocket_center(self):
        """计算口袋中心坐标"""
        all_coords = []
        for residue in self.pocket_residues:
            for atom in residue:
                all_coords.append(atom.get_coord())

        return np.mean(all_coords, axis=0) if all_coords else np.zeros(3)

    def extract_all_features(self):
        """提取所有口袋残基特征"""
        features = {
            'residue_type': [],
            'secondary_structure': [],
            'hydrophobicity': [],
            'sidechain_orientation': [],
            'geometry': [],
            'positions': []
        }

        for residue in self.pocket_residues:
            # 1. 残基类型
            resname = residue.get_resname()
            res_idx = AA_TO_INDEX.get(resname, 20)
            features['residue_type'].append(res_idx)

            # 2. 二级结构
            ss_type = self._get_secondary_structure(residue)
            features['secondary_structure'].append(ss_type)

            # 3. 疏水性
            hydro_score = HYDROPHOBICITY_SCALE.get(resname, 0.0)
            features['hydrophobicity'].append(hydro_score)

            # 4. 侧链朝向
            sidechain_vec = self._compute_sidechain_orientation(residue)
            features['sidechain_orientation'].append(sidechain_vec)

            # 5. 几何特征
            geom_feat = self._compute_geometry_features(residue)
            features['geometry'].append(geom_feat)

            # 6. 残基位置 (Cα原子)
            ca_coord = self._get_ca_coord(residue)
            features['positions'].append(ca_coord)

        # 转换为numpy数组
        for key in features:
            features[key] = np.array(features[key])

        return features

    def _get_secondary_structure(self, residue):
        """获取残基的二级结构"""
        if self.dssp_data is None:
            return 2  # 默认为无规卷曲

        chain_id = residue.parent.id
        res_id = residue.id[1]

        try:
            dssp_key = (chain_id, (' ', res_id, ' '))
            ss = self.dssp_data[dssp_key][2]
            return SS_TO_INDEX.get(ss, 2)
        except:
            return 2

    def _compute_sidechain_orientation(self, residue):
        """计算侧链朝向 (相对于口袋中心)"""
        # 获取Cα和Cβ原子
        ca_coord = self._get_ca_coord(residue)

        # 对于GLY没有Cβ，使用N原子代替
        if 'CB' in residue:
            cb_coord = residue['CB'].get_coord()
        elif 'N' in residue:
            cb_coord = residue['N'].get_coord()
        else:
            return np.array([0, 0, 0])

        # 侧链方向向量
        sidechain_vec = cb_coord - ca_coord

        # 归一化
        norm = np.linalg.norm(sidechain_vec)
        if norm > 0:
            sidechain_vec = sidechain_vec / norm

        return sidechain_vec

    def _compute_geometry_features(self, residue):
        """计算几何特征 (法向量和曲率)"""
        # 获取残基周围的残基
        ca_coord = self._get_ca_coord(residue)
        neighbor_coords = []

        for other_res in self.pocket_residues:
            if other_res == residue:
                continue

            other_ca = self._get_ca_coord(other_res)
            dist = np.linalg.norm(ca_coord - other_ca)

            if dist < 8.0:  # 8Å邻域
                neighbor_coords.append(other_ca)

        if len(neighbor_coords) < 3:
            # 不足3个邻居，无法计算法向量
            return np.concatenate([np.zeros(3), np.zeros(3)])

        # 使用PCA计算局部平面法向量
        neighbor_coords = np.array(neighbor_coords)
        centered = neighbor_coords - ca_coord

        # 协方差矩阵
        cov = np.cov(centered.T)

        # 特征值和特征向量
        eigenvalues, eigenvectors = np.linalg.eig(cov)

        # 最小特征值对应的特征向量是法向量
        idx = np.argmin(eigenvalues)
        normal_vec = eigenvectors[:, idx]

        # 曲率特征 (用特征值表示)
        sorted_eigenvalues = np.sort(eigenvalues)[::-1]
        curvature = sorted_eigenvalues / (np.sum(sorted_eigenvalues) + 1e-6)

        return np.concatenate([normal_vec, curvature])

    def _get_ca_coord(self, residue):
        """获取Cα原子坐标"""
        if 'CA' in residue:
            return residue['CA'].get_coord()
        else:
            # 如果没有CA，使用第一个原子
            return list(residue.get_atoms())[0].get_coord()

    def save_features(self, output_file):
        """保存特征到文件"""
        features = self.extract_all_features()

        np.savez(output_file,
                residue_type=features['residue_type'],
                secondary_structure=features['secondary_structure'],
                hydrophobicity=features['hydrophobicity'],
                sidechain_orientation=features['sidechain_orientation'],
                geometry=features['geometry'],
                positions=features['positions'])

        print(f"口袋特征已保存至: {output_file}")

def extract_pocket_features(pdb_file, ligand_coords, output_file):
    """提取口袋特征的主函数"""
    extractor = PocketFeatureExtractor(pdb_file, ligand_coords)
    extractor.save_features(output_file)
    return extractor.extract_all_features()

if __name__ == "__main__":
    # 示例用法
    pdb_file = "./data/pdbbind/v2020/1a1e/1a1e_protein.pdb"
    ligand_file = "./data/pdbbind/v2020/1a1e/1a1e_ligand.sdf"

    # 读取配体坐标
    from rdkit import Chem
    mol = Chem.SDMolSupplier(ligand_file)[0]
    ligand_coords = mol.GetConformer().GetPositions()

    # 提取特征
    features = extract_pocket_features(pdb_file, ligand_coords,
                                      "./data/1a1e_pocket_features.npz")

    print("特征维度:")
    for key, val in features.items():
        print(f"  {key}: {val.shape}")
