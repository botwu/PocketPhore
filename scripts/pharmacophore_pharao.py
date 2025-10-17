"""
使用Pharao (开源药效团工具) 生成药效团

安装:
conda install -c conda-forge pharao

或者从源码编译:
git clone https://github.com/silicos-it/pharao.git
"""

import subprocess
import json
import os

def generate_pharmacophore_pharao(ligand_sdf, output_phore):
    """
    使用Pharao生成药效团

    Pharao特点:
    - 完全开源
    - 支持多种文件格式
    - 可自定义特征定义
    """
    cmd = [
        'pharao',
        '--dbType', 'SDFFILE',
        '--dbase', ligand_sdf,
        '--pharmacophore',
        '--outType', 'JSON',
        '--out', output_phore
    ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"Pharao生成药效团: {output_phore}")

        # 解析Pharao输出
        with open(output_phore, 'r') as f:
            pharao_data = json.load(f)

        return convert_pharao_to_phoredata(pharao_data)

    except subprocess.CalledProcessError as e:
        print(f"Pharao运行失败: {e}")
        return None

def convert_pharao_to_phoredata(pharao_data):
    """转换Pharao格式为PhoreGen格式"""
    pharmacophore = []

    # Pharao特征映射
    feature_map = {
        'AROM': 'AR',   # 芳香环
        'HDON': 'HD',   # 氢键供体
        'HACC': 'HA',   # 氢键受体
        'POSC': 'PO',   # 正电荷
        'NEGC': 'NE',   # 负电荷
        'HYBL': 'HY',   # 疏水
        'LIPO': 'HY'    # 亲脂
    }

    for feature in pharao_data.get('pharmacophore', []):
        feat_type = feature['type']
        coords = feature['coordinates']

        if feat_type in feature_map:
            pharmacophore.append({
                'type': feature_map[feat_type],
                'position': coords,
                'radius': feature.get('radius', 1.0)
            })

    return pharmacophore
