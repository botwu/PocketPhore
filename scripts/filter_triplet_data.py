"""
三元组数据质量过滤与验证

过滤标准:
1. 配体分子量 < 800 Da
2. 药效团特征数量 3-15个
3. 口袋残基数量 10-100个
4. 结构完整性检查
5. 配体-口袋距离合理性
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import pickle
from tqdm import tqdm

class TripletQualityFilter:
    def __init__(self, config=None):
        # 默认过滤标准
        self.config = config or {
            'min_mw': 200,           # 最小分子量
            'max_mw': 800,           # 最大分子量
            'min_phore_features': 3, # 最少药效团特征
            'max_phore_features': 15,# 最多药效团特征
            'min_pocket_residues': 10,  # 最少口袋残基
            'max_pocket_residues': 100, # 最多口袋残基
            'max_hba': 10,           # 最多氢键受体
            'max_hbd': 5,            # 最多氢键供体
            'max_rotatable_bonds': 10,  # 最多可旋转键
            'min_ligand_pocket_dist': 2.0,  # 最小配体-口袋距离
            'max_ligand_pocket_dist': 15.0  # 最大配体-口袋距离
        }

    def check_ligand_quality(self, ligand_mol):
        """检查配体质量"""
        if ligand_mol is None:
            return False, "配体分子为None"

        # 分子量
        mw = Descriptors.MolWt(ligand_mol)
        if not (self.config['min_mw'] <= mw <= self.config['max_mw']):
            return False, f"分子量超出范围: {mw:.1f}"

        # Lipinski规则
        hba = Lipinski.NumHAcceptors(ligand_mol)
        hbd = Lipinski.NumHDonors(ligand_mol)
        rot_bonds = Lipinski.NumRotatableBonds(ligand_mol)

        if hba > self.config['max_hba']:
            return False, f"氢键受体过多: {hba}"

        if hbd > self.config['max_hbd']:
            return False, f"氢键供体过多: {hbd}"

        if rot_bonds > self.config['max_rotatable_bonds']:
            return False, f"可旋转键过多: {rot_bonds}"

        # 检查是否有3D坐标
        try:
            conf = ligand_mol.GetConformer()
            if conf.GetNumAtoms() == 0:
                return False, "无3D坐标"
        except:
            return False, "无法获取构象"

        return True, "通过"

    def check_pharmacophore_quality(self, pharmacophore):
        """检查药效团质量"""
        if not pharmacophore:
            return False, "药效团为空"

        n_features = len(pharmacophore)

        if not (self.config['min_phore_features'] <= n_features <= self.config['max_phore_features']):
            return False, f"药效团特征数量异常: {n_features}"

        # 检查特征类型的多样性
        feature_types = set([f['type'] for f in pharmacophore])
        if len(feature_types) < 2:
            return False, "药效团特征类型单一"

        # 检查坐标有效性
        for feat in pharmacophore:
            pos = feat['position']
            if len(pos) != 3:
                return False, "药效团坐标维度错误"

            if np.any(np.isnan(pos)) or np.any(np.isinf(pos)):
                return False, "药效团坐标包含NaN或Inf"

        return True, "通过"

    def check_pocket_quality(self, pocket_features):
        """检查口袋质量"""
        if pocket_features is None:
            return False, "口袋特征为None"

        n_residues = len(pocket_features.get('residue_type', []))

        if not (self.config['min_pocket_residues'] <= n_residues <= self.config['max_pocket_residues']):
            return False, f"口袋残基数量异常: {n_residues}"

        # 检查特征完整性
        required_keys = ['residue_type', 'secondary_structure', 'hydrophobicity',
                        'sidechain_orientation', 'geometry', 'positions']

        for key in required_keys:
            if key not in pocket_features:
                return False, f"缺少特征: {key}"

            if len(pocket_features[key]) != n_residues:
                return False, f"特征维度不匹配: {key}"

        return True, "通过"

    def check_spatial_relationship(self, ligand_coords, pocket_positions):
        """检查配体-口袋空间关系"""
        # 计算配体中心
        ligand_center = np.mean(ligand_coords, axis=0)

        # 计算到口袋残基的距离
        distances = np.linalg.norm(pocket_positions - ligand_center, axis=1)

        min_dist = np.min(distances)
        mean_dist = np.mean(distances)

        if min_dist < self.config['min_ligand_pocket_dist']:
            return False, f"配体与口袋过近: {min_dist:.2f}Å"

        if mean_dist > self.config['max_ligand_pocket_dist']:
            return False, f"配体与口袋过远: {mean_dist:.2f}Å"

        return True, "通过"

    def filter_triplet(self, triplet):
        """过滤单个三元组"""
        checks = []

        # 1. 检查配体
        ligand_mol = triplet['ligand']['mol']
        ligand_ok, ligand_msg = self.check_ligand_quality(ligand_mol)
        checks.append(('ligand', ligand_ok, ligand_msg))

        # 2. 检查药效团
        pharmacophore = triplet['pharmacophore']
        phore_ok, phore_msg = self.check_pharmacophore_quality(pharmacophore)
        checks.append(('pharmacophore', phore_ok, phore_msg))

        # 3. 检查口袋 (如果有)
        if 'pocket_features' in triplet:
            pocket_ok, pocket_msg = self.check_pocket_quality(triplet['pocket_features'])
            checks.append(('pocket', pocket_ok, pocket_msg))

            # 4. 检查空间关系
            if pocket_ok and ligand_ok:
                ligand_coords = np.array(triplet['ligand']['coords'])
                pocket_positions = triplet['pocket_features']['positions']

                spatial_ok, spatial_msg = self.check_spatial_relationship(
                    ligand_coords, pocket_positions)
                checks.append(('spatial', spatial_ok, spatial_msg))

        # 汇总结果
        all_passed = all([check[1] for check in checks])

        return all_passed, checks

    def filter_dataset(self, triplets, verbose=True):
        """批量过滤数据集"""
        filtered_triplets = []
        filter_stats = {
            'total': len(triplets),
            'passed': 0,
            'failed_reasons': {}
        }

        for triplet in tqdm(triplets, desc="过滤三元组"):
            passed, checks = self.filter_triplet(triplet)

            if passed:
                filtered_triplets.append(triplet)
                filter_stats['passed'] += 1
            else:
                # 记录失败原因
                for check_name, check_ok, check_msg in checks:
                    if not check_ok:
                        reason = f"{check_name}: {check_msg}"
                        filter_stats['failed_reasons'][reason] = \
                            filter_stats['failed_reasons'].get(reason, 0) + 1

        if verbose:
            self.print_filter_stats(filter_stats)

        return filtered_triplets, filter_stats

    def print_filter_stats(self, stats):
        """打印过滤统计信息"""
        print("\n" + "="*60)
        print("数据过滤统计")
        print("="*60)
        print(f"总数据量: {stats['total']}")
        print(f"通过过滤: {stats['passed']} ({stats['passed']/stats['total']*100:.1f}%)")
        print(f"未通过: {stats['total'] - stats['passed']}")
        print("\n失败原因分布:")

        # 按失败次数排序
        sorted_reasons = sorted(stats['failed_reasons'].items(),
                               key=lambda x: x[1], reverse=True)

        for reason, count in sorted_reasons[:10]:  # 显示前10个原因
            print(f"  {reason}: {count}")

        print("="*60 + "\n")

def validate_and_save_dataset(input_file, output_file, config=None):
    """验证并保存过滤后的数据集"""
    # 加载原始数据
    print(f"加载数据: {input_file}")
    with open(input_file, 'rb') as f:
        triplets = pickle.load(f)

    # 创建过滤器
    filter = TripletQualityFilter(config)

    # 过滤数据
    filtered_triplets, stats = filter.filter_dataset(triplets)

    # 保存过滤后的数据
    with open(output_file, 'wb') as f:
        pickle.dump(filtered_triplets, f)

    print(f"过滤后数据已保存至: {output_file}")

    # 保存统计信息
    stats_file = output_file.replace('.pkl', '_stats.json')
    import json
    with open(stats_file, 'w') as f:
        # 转换为可序列化的格式
        serializable_stats = {
            'total': stats['total'],
            'passed': stats['passed'],
            'pass_rate': stats['passed'] / stats['total'],
            'failed_reasons': stats['failed_reasons']
        }
        json.dump(serializable_stats, f, indent=2)

    return filtered_triplets, stats

if __name__ == "__main__":
    # 示例：过滤PDBbind数据
    input_file = './data/triplets.pkl'
    output_file = './data/triplets_filtered.pkl'

    # 自定义过滤配置
    custom_config = {
        'min_mw': 200,
        'max_mw': 800,
        'min_phore_features': 3,
        'max_phore_features': 15,
        'min_pocket_residues': 10,
        'max_pocket_residues': 100,
        'max_hba': 10,
        'max_hbd': 5,
        'max_rotatable_bonds': 10,
        'min_ligand_pocket_dist': 2.0,
        'max_ligand_pocket_dist': 15.0
    }

    filtered_triplets, stats = validate_and_save_dataset(
        input_file, output_file, config=custom_config
    )
