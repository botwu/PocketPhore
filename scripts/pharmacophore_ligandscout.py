"""
使用LigandScout命令行工具生成药效团

需要安装: LigandScout (商业软件，有学术免费版)
下载: https://www.inteligand.com/ligandscout/
"""

import subprocess
import os
import xml.etree.ElementTree as ET

def generate_pharmacophore_ligandscout(pdb_file, ligand_file, output_phore):
    """
    使用LigandScout生成药效团

    LigandScout特点：
    - 自动识别11种药效团特征
    - 考虑蛋白-配体相互作用
    - 输出标准.pml格式
    """
    cmd = [
        'idbgen',  # LigandScout命令行工具
        '-p', pdb_file,
        '-l', ligand_file,
        '-o', output_phore,
        '--pharmacophore-type', 'structure-based',  # 基于结构的药效团
        '--include-exclusion-volumes',  # 包含排斥体积
        '--feature-tolerance', '1.0'  # 特征容忍度
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True)
        print(f"成功生成药效团: {output_phore}")
        return parse_ligandscout_pharmacophore(output_phore)
    except subprocess.CalledProcessError as e:
        print(f"LigandScout运行失败: {e}")
        return None

def parse_ligandscout_pharmacophore(pml_file):
    """解析LigandScout .pml文件"""
    tree = ET.parse(pml_file)
    root = tree.getroot()

    pharmacophore = []

    # LigandScout特征类型映射
    feature_map = {
        'HBD': 'HD',  # 氢键供体
        'HBA': 'HA',  # 氢键受体
        'AR': 'AR',   # 芳香环
        'PI': 'PO',   # 正电离子
        'NI': 'NE',   # 负电离子
        'H': 'HY',    # 疏水
        'XI': 'XV'    # 排斥体积
    }

    for feature in root.findall('.//pharmacophore/point'):
        feat_type = feature.get('name')
        position = feature.find('position')

        x = float(position.get('x3'))
        y = float(position.get('y3'))
        z = float(position.get('z3'))

        tolerance = feature.find('tolerance')
        radius = float(tolerance.get('radius')) if tolerance is not None else 1.0

        if feat_type in feature_map:
            pharmacophore.append({
                'type': feature_map[feat_type],
                'position': [x, y, z],
                'radius': radius
            })

    return pharmacophore
