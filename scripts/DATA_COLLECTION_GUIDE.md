# 蛋白-配体-药效团三元组数据收集指南

## 目录
1. [数据源](#数据源)
2. [环境准备](#环境准备)
3. [快速开始](#快速开始)
4. [详细步骤](#详细步骤)
5. [数据格式](#数据格式)
6. [FAQ](#faq)

---

## 数据源

### 推荐数据集

| 数据库 | 数量 | 下载链接 | 优先级 |
|--------|------|---------|--------|
| **PDBbind** | ~6,000 | http://www.pdbbind.org.cn/download.php | ⭐⭐⭐ |
| **CrossDocked2020** | ~22M | https://github.com/gnina/models | ⭐⭐⭐ |
| **Binding MOAD** | ~35,000 | http://bindingmoad.org | ⭐⭐ |

### PDBbind目录结构
```
pdbbind/
├── v2020/
│   ├── index/
│   │   └── INDEX_refined_data.2020
│   ├── 1a1e/
│   │   ├── 1a1e_protein.pdb
│   │   ├── 1a1e_ligand.sdf
│   │   └── 1a1e_pocket.pdb
│   ├── 1a28/
│   └── ...
```

---

## 环境准备

### 1. 基础依赖

```bash
# 创建conda环境
conda create -n pocketphore python=3.9
conda activate pocketphore

# 安装核心包
pip install rdkit-pypi biopython scipy numpy tqdm

# 安装DSSP (二级结构计算)
conda install -c salilab dssp
```

### 2. 可选工具 (提升药效团质量)

#### LigandScout (商业软件，有学术免费版)
```bash
# 下载地址: https://www.inteligand.com/ligandscout/
# 安装后添加到PATH
export PATH=$PATH:/path/to/LigandScout/bin
```

#### Pharao (开源)
```bash
conda install -c conda-forge pharao
```

---

## 快速开始

### 方法1: 一键运行完整Pipeline (推荐)

```bash
# 下载PDBbind数据集
# http://www.pdbbind.org.cn/download.php

# 运行完整pipeline
python scripts/full_pipeline.py \
    --pdbbind_dir ./data/pdbbind \
    --output_dir ./data/processed \
    --n_workers 8

# 输出文件:
# - triplets_raw.pkl: 原始三元组
# - triplets_filtered.pkl: 过滤后的三元组
# - train.pkl, val.pkl, test.pkl: 划分后的数据集
```

### 方法2: 分步执行

```bash
# Step 1: 收集原始数据
python scripts/collect_triplet_data.py

# Step 2: 质量过滤
python scripts/filter_triplet_data.py \
    --input ./data/triplets_raw.pkl \
    --output ./data/triplets_filtered.pkl

# Step 3: 数据集划分 (在Python中)
python -c "
import pickle
import numpy as np

with open('./data/triplets_filtered.pkl', 'rb') as f:
    data = pickle.load(f)

n = len(data)
idx = np.random.permutation(n)
splits = {
    'train': [data[i] for i in idx[:int(n*0.8)]],
    'val': [data[i] for i in idx[int(n*0.8):int(n*0.9)]],
    'test': [data[i] for i in idx[int(n*0.9):]]
}

for name, split in splits.items():
    with open(f'./data/{name}.pkl', 'wb') as f:
        pickle.dump(split, f)
    print(f'{name}: {len(split)} samples')
"
```

---

## 详细步骤

### Step 1: 蛋白-配体复合物提取

```python
from rdkit import Chem
from Bio.PDB import PDBParser

# 读取配体
ligand_mol = Chem.SDMolSupplier('1a1e_ligand.sdf')[0]
ligand_coords = ligand_mol.GetConformer().GetPositions()

# 读取蛋白
parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', '1a1e_protein.pdb')
```

### Step 2: 口袋残基提取

```python
from scripts.extract_pocket_features import PocketFeatureExtractor

# 提取配体10Å周围的残基
extractor = PocketFeatureExtractor(
    '1a1e_protein.pdb',
    ligand_coords,
    pocket_cutoff=10.0
)

# 提取特征
pocket_features = extractor.extract_all_features()

# 特征包括:
# - residue_type: (N_res,) 残基类型索引
# - secondary_structure: (N_res,) 二级结构类型
# - hydrophobicity: (N_res,) 疏水性评分
# - sidechain_orientation: (N_res, 3) 侧链朝向
# - geometry: (N_res, 6) 几何特征 (法向量+曲率)
# - positions: (N_res, 3) 残基Cα坐标
```

### Step 3: 药效团生成

#### 方法A: 使用RDKit (默认)
```python
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig

fdef_name = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdef_name)
features = factory.GetFeaturesForMol(ligand_mol)

pharmacophore = []
for feat in features:
    pharmacophore.append({
        'type': feat.GetFamily(),  # HD, HA, AR, HY, PO, NE
        'position': [feat.GetPos().x, feat.GetPos().y, feat.GetPos().z],
        'radius': 1.0
    })
```

#### 方法B: 使用LigandScout (推荐)
```bash
idbgen -p 1a1e_protein.pdb \
       -l 1a1e_ligand.sdf \
       -o 1a1e.pml \
       --pharmacophore-type structure-based
```

### Step 4: 数据质量过滤

```python
from scripts.filter_triplet_data import TripletQualityFilter

# 设置过滤标准
filter = TripletQualityFilter(config={
    'min_mw': 200,              # 最小分子量
    'max_mw': 800,              # 最大分子量
    'min_phore_features': 3,    # 最少药效团特征
    'max_phore_features': 15,   # 最多药效团特征
    'min_pocket_residues': 10,  # 最少口袋残基
    'max_pocket_residues': 100, # 最多口袋残基
})

# 过滤数据
filtered_data, stats = filter.filter_dataset(triplets)

# 输出统计信息:
# 总数据量: 6000
# 通过过滤: 4850 (80.8%)
# 未通过: 1150
```

---

## 数据格式

### 三元组数据结构

```python
triplet = {
    'pdb_id': '1a1e',

    'ligand': {
        'mol': <RDKit Mol对象>,
        'smiles': 'CCO...',
        'coords': [[x1,y1,z1], [x2,y2,z2], ...],  # (N_atoms, 3)
        'file': '1a1e_ligand.sdf'
    },

    'protein': {
        'pocket_file': '1a1e_pocket.pdb',
        'full_file': '1a1e_protein.pdb'
    },

    'pharmacophore': [
        {'type': 'HD', 'position': [x, y, z], 'radius': 1.0},
        {'type': 'HA', 'position': [x, y, z], 'radius': 1.0},
        ...
    ],

    'pocket_features': {
        'residue_type': np.array([0, 5, 12, ...]),  # (N_res,)
        'secondary_structure': np.array([0, 0, 1, ...]),  # (N_res,)
        'hydrophobicity': np.array([1.8, -0.4, ...]),  # (N_res,)
        'sidechain_orientation': np.array([[...], ...]),  # (N_res, 3)
        'geometry': np.array([[...], ...]),  # (N_res, 6)
        'positions': np.array([[...], ...])  # (N_res, 3)
    }
}
```

### 药效团特征类型映射

| PhoreGen类型 | 含义 | RDKit对应 | LigandScout对应 |
|-------------|------|-----------|----------------|
| HD | 氢键供体 | Donor | HBD |
| HA | 氢键受体 | Acceptor | HBA |
| AR | 芳香环 | Aromatic | AR |
| HY | 疏水 | Hydrophobe | H |
| PO | 正电 | PosIonizable | PI |
| NE | 负电 | NegIonizable | NI |
| XV | 排斥体积 | - | XI |

---

## FAQ

### Q1: PDBbind数据集下载很慢怎么办?
**A**: 使用国内镜像或代理，或者使用CrossDocked2020数据集替代。

### Q2: DSSP安装失败怎么办?
**A**:
```bash
# macOS
brew install dssp

# Linux
sudo apt-get install dssp

# 或使用conda
conda install -c salilab dssp
```

### Q3: 没有LigandScout license怎么办?
**A**: 使用RDKit生成药效团（默认方式），质量略低但免费。或申请LigandScout学术免费版。

### Q4: 内存不足怎么办?
**A**:
- 减少`n_workers`并行数
- 分批处理数据
- 使用交换空间

```bash
# 分批处理示例
python full_pipeline.py --pdbbind_dir ./data/pdbbind \
    --output_dir ./data/batch1 \
    --pdb_ids_file ./pdb_ids_batch1.txt
```

### Q5: 如何验证数据质量?
**A**:
```python
import pickle
from rdkit import Chem

# 加载数据
with open('./data/train.pkl', 'rb') as f:
    train_data = pickle.load(f)

# 检查第一个样本
sample = train_data[0]
print(f"PDB ID: {sample['pdb_id']}")
print(f"配体原子数: {sample['ligand']['mol'].GetNumAtoms()}")
print(f"药效团特征数: {len(sample['pharmacophore'])}")
print(f"口袋残基数: {len(sample['pocket_features']['residue_type'])}")

# 可视化 (使用py3Dmol)
import py3Dmol
view = py3Dmol.view(width=800, height=600)
view.addModel(open(sample['protein']['pocket_file']).read(), 'pdb')
view.setStyle({'cartoon': {'color': 'spectrum'}})
view.zoomTo()
view.show()
```

### Q6: 数据集推荐使用规模?
**A**:
- **最小训练集**: 1,000个三元组 (概念验证)
- **推荐训练集**: 5,000-10,000个三元组 (较好性能)
- **理想训练集**: 20,000+个三元组 (最佳性能)

### Q7: 如何加速数据收集?
**A**:
```bash
# 使用更多CPU核心
python full_pipeline.py --n_workers 16

# 跳过DSSP计算 (如果不需要二级结构)
python full_pipeline.py --skip_dssp

# 使用GPU加速特征提取 (如果有)
python full_pipeline.py --use_gpu
```

---

## 预期数据统计

基于PDBbind v2020精炼集：

| 统计项 | 数值 |
|--------|------|
| 总复合物数 | ~6,000 |
| 过滤后数量 | ~4,800 (80%) |
| 平均配体原子数 | 25-30 |
| 平均药效团特征数 | 6-8 |
| 平均口袋残基数 | 35-45 |
| 磁盘占用 | ~10GB (原始) + ~5GB (处理后) |
| 处理时间 | ~2-4小时 (8核CPU) |

---

## 引用

如果使用PDBbind数据集，请引用：
```
Wang, R., Fang, X., Lu, Y., & Wang, S. (2004).
The PDBbind database: collection of binding affinities for protein− ligand complexes with known three-dimensional structures.
Journal of medicinal chemistry, 47(12), 2977-2980.
```

---

## 联系方式

遇到问题请提Issue或联系：
- GitHub: https://github.com/yourname/PhoreGen
- Email: your.email@example.com
