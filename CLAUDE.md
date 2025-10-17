# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PhoreGen is a pharmacophore-oriented 3D molecular generation framework that generates molecules aligned with pharmacophore models using diffusion-based deep learning. It employs asynchronous perturbations on both atomic and bond information with message-passing mechanisms that incorporate ligand-pharmacophore mapping knowledge.

## Environment Setup

### Standard Installation (Linux with CUDA)
```bash
conda env create -f phoregen_env.yml
conda activate phoregen
```

### macOS/CPU-only Installation
The conda environment file is designed for Linux with CUDA. For macOS or CPU-only systems:
```bash
conda create -n phoregen python=3.9
conda activate phoregen
pip install torch torchvision torchaudio
pip install rdkit-pypi torch-geometric pyyaml easydict tensorboardx wandb scipy matplotlib tqdm pandas lmdb
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-2.2.0+cpu.html
conda install -c conda-forge openbabel
```

**Key Dependencies**: Python 3.9, PyTorch 1.12.1+, PyTorch Geometric 2.1.0+, RDKit 2022.9.5, OpenBabel 3.1.1

## Core Commands

### Molecular Generation (Sampling)
Generate molecules from pharmacophore models:
```bash
python sample_all.py --num_samples 100 \
    --outdir ./results/test \
    --phore_file_list ./data/phores_for_sampling/file_index.json \
    --check_point ./ckpt/zinc_trained.pt \
    --batch_size 30
```

Key sampling arguments:
- `--num_samples`: Number of molecules per pharmacophore (default: 1000)
- `--batch_size`: Samples per generation step (default: 30)
- `--check_point`: Model checkpoint path (required)
- `--config`: Config file (default: `train_dock-cpx-phore.yml`)
- `--add_edge`: Edge prediction mode ('predicted' or 'fully')
- `--pos_guidance_opt`: JSON string for position guidance (e.g., atom proximity constraints)
- `--sample_nodes_mode`: 'uniform' or 'normal' for node count sampling
- `--normal_scale`: Scale parameter for normal distribution sampling

Batch sampling with shell script:
```bash
bash sample.sh
```

### Training

Pre-training on LigPhore dataset:
```bash
python train.py --config ./configs/train_lig-phore.yml
```

Fine-tuning on CpxPhore and DockPhore:
```bash
python train.py --config ./configs/train_dock-cpx-phore.yml
```

## Architecture Overview

### Diffusion Process

PhoreGen uses a three-component asynchronous diffusion process:
1. **Position Diffusion**: Continuous Gaussian noise on 3D coordinates
2. **Atom Type Diffusion**: Discrete categorical diffusion with "mask" state
3. **Bond Type Diffusion**: Discrete categorical diffusion with "absorb" state

The diffusion process runs for 1000 timesteps with segment-based scheduling for bonds (600 steps, then 400 steps with different noise scales).

### Core Components

**PhoreDiff** (`models/diffusion.py`)
- Main diffusion model orchestrating the generation process
- Handles both atomic and bond information with asynchronous perturbations
- Manages three separate transition models for position, atoms, and bonds
- Loss composition: position (weight=1) + node type (weight=100) + edge type (weight=100)

**Denoiser Network** (`models/uni_denoiser.py`)
- Graph transformer architecture with 6 layers, 16 attention heads
- KNN-based graph construction (k=32 neighbors)
- Processes ligand atoms and pharmacophore features jointly
- Returns denoised positions, node embeddings, and bond embeddings

**Transition Models** (`models/transition.py`)
- `ContigousTransition`: For continuous position diffusion
- `CategoricalTransition`: Base class for discrete atom/bond diffusion
- `GeneralCategoricalTransition`: Extended categorical transition with multiple initialization strategies

**Data Processing** (`datasets/`)
- `phoregen.py`: Main dataset classes for ligand-pharmacophore pairs
- `get_phore_data.py`: Loads pharmacophore models from .phore files
- `generate_phorefp.py`: Generates pharmacophore fingerprints
- Supports three datasets: LigPhore, CpxPhore, DockPhore (available at Zenodo)

### Data Flow

1. **Loading**: Pharmacophore `.phore` files + molecular data → `PhoreData`
2. **Preprocessing**: Transform adds noise, creates fully-connected edges
3. **Forward Pass**: Diffusion model adds noise → Denoiser predicts clean states
4. **Sampling**: Reverse diffusion from noise → Molecule reconstruction
5. **Output**: Generated molecules saved as `.sdf` files + SMILES strings

### Configuration Architecture

Configs are hierarchical YAML files with these main sections:
- `model`: Diffusion parameters, network architecture, feature dimensions
- `train`: Training hyperparameters, device, optimizer, scheduler
- `dataset`: Data paths, preprocessing options, dataset name
- `logger`: Output paths, run names, checkpointing

**Critical Configuration Details**:
- `phore_feat_dim`: Automatically adjusted based on dataset:
  - For `zinc_300` and `pdbbind`: +2 features added programmatically
  - Base value: 16 or 18 depending on dataset
- `device`: Must be set to 'cpu' for non-CUDA systems (defaults to 'cuda')
- `data_name`: Controls feature augmentation logic; must match checkpoint training

## Molecular Representation

**Atom Types** (12 classes): B, C, N, O, F, Si, P, S, Cl, Br, I + masked_atom
**Bond Types** (6 classes): 0 (no bond), 1 (single), 2 (double), 3 (triple), 4 (aromatic) + masked_bond
**Pharmacophore Features** (11 types): MB, HD, AR, PO, HA, HY, NE, CV, CR, XB, EX
- MB: Metal Binding
- HD: Hydrogen Bond Donor
- AR: Aromatic Ring
- PO: Positive Ionizable
- HA: Hydrogen Bond Acceptor
- HY: Hydrophobic
- NE: Negative Ionizable
- CV: Covalent
- CR: Chelating Ring
- XB: Halogen Bond
- EX: Exclusion Volume

## Pretrained Models

Download from [Zenodo](https://zenodo.org/records/14404575):
- `zinc_trained.pt`: Pre-trained on LigPhore (ZINC-derived ligands)
- `crossdocked_pdbbind_trained.pt`: Fine-tuned on CpxPhore + DockPhore

Place in `./ckpt/` directory.

## Common Issues and Solutions

### Device Mismatch
**Error**: `AssertionError: Torch not compiled with CUDA enabled`
**Solution**: Set `device: cpu` in config files (`train.device` section)

### Model Parameter Mismatch
**Error**: `size mismatch for phore_embedding.weight`
**Solution**: The `phore_feat_dim` is automatically adjusted for certain datasets. Ensure:
- Config `data_name` matches checkpoint training dataset
- For `zinc_300`/`pdbbind`: config should have `phore_feat_dim: 16` (becomes 18 after +2)
- For other datasets: use `phore_feat_dim: 16` without adjustment

### Deprecated PyTorch Functions
The code contains `torch.cross()` without `dim` argument. This has been fixed but may appear in warnings. Use `torch.cross(pos_ji, pos_ki, dim=-1)` format.

### Slow Sampling Performance
Diffusion sampling is computationally intensive (~4+ minutes per molecule on CPU):
- Use GPU when possible
- Reduce `num_samples` for testing
- Increase `batch_size` if memory allows (GPU only)

## Output Structure

Generated molecules are saved in the output directory:
```
results/test/
├── <phore_name>_samples_all.pt          # PyTorch tensor with all sample data
├── <phore_name>_SMILES_all.txt          # SMILES strings (one per line)
├── time_chain.txt                       # Generation time statistics
└── sdf_results/
    └── <idx>_<phore_name>_<sample>.sdf  # 3D molecular structures
```

## Generating Custom Pharmacophore Models

Use the online tool [AncPhore](https://ancphore.ddtmlab.org/Model) to create pharmacophore models from:
- Protein-ligand complexes
- Individual ligands
- Custom feature points

Export as `.phore` files and add to `phore_file_list` JSON array.

## Code Modification Guidelines

When modifying the diffusion model:
- Loss weights affect training stability: [1, 100, 100] for [pos, node, edge] is well-tuned
- Time step sampling is uniform over [0, 999] - avoid biasing toward start/end
- Pharmacophore features are embedded jointly with ligand atoms - maintain this structure
- The fully-connected edge representation is crucial for bond prediction - don't skip it
