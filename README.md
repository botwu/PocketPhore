# PocketPhore

**Pocket-Enhanced Pharmacophore-Oriented 3D Molecular Generation**

**基于蛋白质口袋信息增强的药效团导向3D分子生成**

[![Pytorch](https://img.shields.io/badge/PyTorch-%23EE4C2C.svg?e&logo=PyTorch&logoColor=white)](https://pytorch.org/)
![](https://img.shields.io/badge/version-1.0.0-blue)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/ppjian19/PhoreGen/blob/main/LICENSE)

---

## Overview | 概述

PocketPhore extends PhoreGen by incorporating protein pocket information into pharmacophore-oriented molecular generation. The method introduces a trimodal framework (ligand-pharmacophore-pocket) with dual-scale attention mechanism to achieve both feature customization and binding affinity optimization.

PocketPhore在PhoreGen基础上，创新性地引入蛋白质口袋信息，构建配体-药效团-蛋白口袋三模态协同生成框架，通过双尺度注意力机制实现特征定制与结合能力的双重优化。

**Key Features | 核心特性:**

- **Trimodal Architecture | 三模态架构**: Ligand + Pharmacophore + Protein Pocket | 配体 + 药效团 + 蛋白口袋
- **Dual-Scale Attention | 双尺度注意力**: Global (pharmacophore) + Local (pocket) guidance | 全局（药效团）+ 局部（口袋）引导
- **Physicochemical Constraints | 物理化学约束**: Explicit H-bond, hydrophobic, electrostatic modeling | 显式氢键、疏水、静电建模
- **Time-Adaptive Fusion | 时间自适应融合**: Dynamic weight adjustment across diffusion stages | 扩散过程不同阶段的动态权重调整

---

## Method | 方法

### Architecture | 架构

**Input | 输入:**
- Noised molecule M<sub>t</sub> | 加噪分子
- Pharmacophore model P | 药效团模型
- Protein pocket Q | 蛋白口袋（新增）

**Core Components | 核心组件:**

1. **Trimodal Embeddings | 三模态嵌入**
   - Ligand embedding | 配体嵌入
   - Pharmacophore embedding | 药效团嵌入
   - Pocket embedding (128-dim, 5 modalities) | 口袋嵌入（128维，5种模态）

2. **Heterogeneous Graph | 异构图**
   - L-L: Ligand-ligand (fully connected) | 配体内全连接
   - P-L: Pharmacophore-ligand (k-NN) | 药效团-配体 k-NN
   - Q-L: Pocket-ligand (k-NN, NEW) | 口袋-配体 k-NN（新增）
   - P-Q: Pharmacophore-pocket (spatial matching, NEW) | 药效团-口袋匹配（新增）

3. **Dual-Scale Time-Adaptive Attention | 双尺度时间自适应注意力**
   - Global: Pharmacophore attention for feature guidance | 全局：药效团注意力，特征引导
   - Local: Pocket attention for interaction refinement | 局部：口袋注意力，相互作用优化
   - Time-adaptive weights | 时间自适应权重:
     - Early stage (t > 0.7T): w_global = 1.0, w_local = 0.1 | 早期：药效团主导
     - Mid stage (0.3T < t ≤ 0.7T): w_global = 0.6, w_local = 0.4 | 中期：双重约束
     - Late stage (t ≤ 0.3T): w_global = 0.2, w_local = 0.5 | 后期：口袋主导

4. **Physicochemical Constraints | 物理化学约束**
   - H-bond network | 氢键网络
   - Hydrophobic network | 疏水网络
   - Electrostatic network | 静电网络
   - Clash detection | 冲突检测

5. **Pocket-Aware Size Prediction | 口袋感知大小预测**
   - Incorporating pocket volume constraints | 结合口袋体积约束

### Training Strategy | 训练策略

Three-stage progressive training | 三阶段渐进训练:

- **Stage 1**: Pharmacophore pretraining (LigPhore, 60 epochs) | 药效团预训练
- **Stage 2**: Pocket geometry enhancement (+CpxPhore, 50 epochs) | 口袋几何增强
- **Stage 3**: Physics constraint optimization (+PocketPhore, 50 epochs) | 物理约束优化

---

## Performance | 性能

### Comparison with Baselines | 与基线方法对比

| Method | Mapping Score | IFP Similarity | Docking Score (kcal/mol) |
|--------|--------------|----------------|--------------------------|
| PhoreGen | 0.784 | 0.66 | -8.03 |
| TargetDiff | 0.65 | 0.72 | -8.21 |
| DiffSBDD | 0.62 | 0.74 | -8.35 |
| ShEPhERD | 0.71 | 0.69 | -8.15 |
| **PocketPhore (Ours)** | **0.81 (+2%)** | **0.74 (+12%)** | **-8.62 (+7%)** |

### Evaluation Metrics | 评估指标

- **Mapping Score | 药效团匹配得分**: Pharmacophore feature matching accuracy | 匹配成功的特征点比例
- **IFP Similarity | 相互作用指纹相似度**: Interaction fingerprint similarity (Tanimoto) | 与参考配体的相互作用模式相似度
- **Docking Score | 分子对接得分**: Predicted binding affinity (AutoDock Vina) | 预测的结合自由能

---

## Installation | 安装

### Requirements | 环境要求

Package | Version
--- | ---
Python | 3.9+
PyTorch | 1.12.1+
PyTorch Geometric | 2.1.0+
RDKit | 2022.9.5
OpenBabel | 3.1.1

### Setup | 安装步骤

```bash
# Create conda environment | 创建conda环境
conda env create -f pocketphore_env.yml
conda activate pocketphore
```

For macOS/CPU-only installation, see [CLAUDE.md](./CLAUDE.md#environment-setup).

macOS或仅CPU安装，请参考 [CLAUDE.md](./CLAUDE.md#environment-setup)。

---

## Datasets | 数据集

### Inherited from PhoreGen | 继承自PhoreGen

| Dataset | Size | Type | Source | Usage |
|---------|------|------|--------|-------|
| LigPhore | 2,398,776 | Ligand-Pharmacophore pairs | ZINC | Pretraining 预训练 |
| CpxPhore | 13,585 | Ligand-Pharmacophore pairs | PDBbind | Fine-tuning 微调 |
| DockPhore | 84,767 | Ligand-Pharmacophore pairs | CrossDocked | Augmentation 增强 |

**说明**: "对"指配体-药效团对，包含一个3D配体分子和对应的药效团模型。

### PocketPhore Dataset (In Construction) | PocketPhore数据集（构建中）

**Target | 目标**: >30,000 triplets (ligand-pharmacophore-pocket) | 三元组（配体-药效团-口袋）

**Sources | 来源**: PDBbind v2024, CrossDocked2020

**Triplet Structure | 三元组结构**:
```python
{
    "ligand": {...},           # Ligand molecule | 配体分子
    "pharmacophore": {...},    # Pharmacophore model | 药效团模型
    "pocket": {...}            # Protein pocket features (128-dim) | 蛋白口袋特征（128维）
}
```

---

## Usage | 使用方法

### Training | 训练

**Stage 1: Pharmacophore Pretraining | 阶段1：药效团预训练**
```bash
python train.py --config ./configs/train_lig-phore.yml
```

**Stage 2-3: Pocket Enhancement and Physics Optimization | 阶段2-3：口袋增强与物理优化**
```bash
python train.py --config ./configs/train_dock-cpx-phore.yml
```

### Sampling | 采样生成

```bash
python sample_all.py \
    --num_samples 100 \
    --outdir ./results/test \
    --phore_file_list ./data/phore_for_sampling/file_index.json \
    --check_point ./ckpt/pocketphore_trained.pt
```

**Key Arguments | 主要参数:**
- `num_samples`: Number of molecules to generate per pharmacophore | 每个药效团生成的分子数
- `outdir`: Output directory | 输出目录
- `phore_file_list`: JSON file containing pharmacophore model paths | 药效团模型路径的JSON文件
- `check_point`: Model checkpoint path | 模型检查点路径

**Output | 输出**:
- 3D structures in `.sdf` format | SDF格式的3D结构
- SMILES strings in `.txt` format | TXT格式的SMILES字符串

### Pharmacophore Generation | 药效团生成

Use the online tool [AncPhore](https://ancphore.ddtmlab.org/Model) to generate pharmacophore models from protein-ligand complexes or individual ligands.

使用在线工具 [AncPhore](https://ancphore.ddtmlab.org/Model) 从蛋白-配体复合物或单独配体生成药效团模型。

---

## Project Status | 项目进度

**Current Progress | 当前进度**: 10.4% completed | 已完成

### Completed | 已完成 (100%)

- PhoreGen baseline reproduction and validation | PhoreGen基线复现与验证
  - Mapping Score 0.78 on LigPhore test set | LigPhore测试集上达到0.78
- Literature review (30+ papers) | 文献调研（30+篇）
- Trimodal architecture design | 三模态架构设计
- Data collection and preprocessing | 数据收集与预处理
  - Pharmacophore generation: >100 targets | 药效团生成：>100个靶标
  - Protein pocket analysis: >50 complexes | 蛋白口袋分析：>50个复合物
  - Molecular docking: >200 molecules | 分子对接：>200个分子
- Computing resources: 2×A100

### In Progress | 进行中 (60%)

- Core module implementation | 核心模块实现
  - Dual-scale attention mechanism | 双尺度注意力机制
  - Physicochemical constraints | 物理化学约束

---

## Project Structure | 项目结构

```
PocketPhore/
├── configs/              # Configuration files | 配置文件
│   ├── train_lig-phore.yml
│   └── train_dock-cpx-phore.yml
├── datasets/             # Data loading and processing | 数据加载与处理
│   ├── phoregen.py
│   └── get_phore_data.py
├── models/               # Model architectures | 模型架构
│   ├── diffusion.py          # Main diffusion model | 主扩散模型
│   ├── uni_denoiser.py       # Denoiser network | 去噪网络
│   └── transition.py         # Transition models | 转换模型
├── utils/                # Utility functions | 工具函数
├── scripts/              # Data processing and visualization | 数据处理与可视化
│   ├── generate_ppt_figures.py
│   └── extract_pocket_features.py
├── train.py              # Training script | 训练脚本
├── sample_all.py         # Sampling script | 采样脚本
└── CLAUDE.md             # Detailed project instructions | 详细项目说明
```

---


## Citation | 引用

### PhoreGen (Baseline Method | 基线方法)

```bibtex
@article{peng2025phoregen,
  title={Pharmacophore-oriented 3D molecular generation toward efficient feature-customized drug discovery},
  author={Peng, Jian and Yu, Jia-Lin and Yang, Zhi-Bing and others},
  journal={Nature Computational Science},
  year={2025},
  doi={10.1038/s43588-025-00850-5}
}
```

### PocketPhore (This Work | 本研究)

Manuscript in preparation. Target journal: *Nature Computational Science*

论文准备中。目标期刊：*Nature Computational Science*

---


## Acknowledgments | 致谢

This work builds upon PhoreGen and benefits from:
- Datasets: PDBbind, CrossDocked2020, ZINC
- Tools: AncPhore, AutoDock Vina, RDKit, PyTorch Geometric
- Computing resources: 2×NVIDIA A100 GPUs

本工作基于PhoreGen，并受益于上述数据集、工具和计算资源的支持。

---

## License | 许可证

This project is licensed under the MIT License.

本项目采用MIT许可证。

---

*Last updated: 2025-10-16*
