# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the EC-CCM (Error Control - Chemistry Coordinate Mapping) release project containing mature CCM implementations for multiple OpenFOAM versions. The project implements custom chemistry models that group CFD cells with similar thermochemical states (Y, T, P) to share reaction rates, with optimization for computational efficiency.

## Directory Structure

```
EC-CCM/
├── OpenFOAM-7/
│   ├── code/                # OpenFOAM-7 CCM implementation
│   └── Sandia/              # OpenFOAM-7 test case
├── OpenFOAM-10/
│   ├── code/                # OpenFOAM-10 CCM implementation
│   └── Sandia/              # OpenFOAM-10 test case
└── OpenFOAM-13/
    ├── code/                # OpenFOAM-13 CCM implementation
    └── Sandia/              # OpenFOAM-13 test case
```

## Core Components

### Code Structure by OpenFOAM Version

**OpenFOAM-7**: Simplified structure with integrated chemistry model
- **chemistryModel/**: Main CCM chemistry model with `makeChemistryModel.H`
- **interface/**: Basic chemistry models and solver types
- **helpers/**: Utility classes for debugging, communication, and combustion calculations

**OpenFOAM-10 & 13**: Extended structure with modular components
- **chemistryModel/**: Main CCM model with reduction/tabulation support
- **odeChemistryModel/**: ODE-based chemistry solver integration
- **interface/**: Chemistry solvers and ODE solver interfaces  
- **helpers/**: Complete utility suite including parallel communication
- **ReactionEntry/** & **RoundRobin/**: Reaction data and communication management

## Build & Test Workflow

### OpenFOAM Version Selection
Choose your OpenFOAM version and activate the corresponding environment:

```bash
# For OpenFOAM-7
source /home/zhouyuchen/OpenFOAM/OpenFOAM-7/OpenFOAM-7/etc/bashrc
cd OpenFOAM-7/code/ && wclean && wmake -j

# For OpenFOAM-10  
source /home/zhouyuchen/OpenFOAM/OpenFOAM-10/OpenFOAM-10/etc/bashrc
cd OpenFOAM-10/code/ && wclean && wmake -j

# For OpenFOAM-13
source /home/zhouyuchen/OpenFOAM/OpenFOAM-13/OpenFOAM-13/etc/bashrc
cd OpenFOAM-13/code/ && wclean && wmake -j
```

### Run Test Cases
```bash
# Navigate to the version-relevant test case directory
cd OpenFOAM-7/Sandia/    # for OpenFOAM-7
cd OpenFOAM-10/Sandia/   # for OpenFOAM-10  
cd OpenFOAM-13/Sandia/   # for OpenFOAM-13

# Then the procedure is similar for all versions:

# Parallel execution (6 cores)
decomposePar -force -fileHandler collated >& log.decomposePar
mpirun -n 6 reactingFoam -parallel -fileHandler collated >& log.reactingFoam&

# Single core execution (use timeout for testing - run less than 1 min)
timeout 50s reactingFoam >& log.reactingFoam&
```

### Configuration
- `system/controlDict` - Loads CCM library
- `constant/chemistryProperties` - Set `method CCM`, configure communication optimization and timing
  - `optimizedCommunication true/false` - Targeted vs broadcast communication
  - `debugTime true/false` - Enable detailed timing analysis

## Version-Specific Features

### OpenFOAM-7 Implementation
- **Status**: ✅ Mature and stable with full support
- **Features**: Complete CCM functionality with ecVars system
- **Testing**: On-going validation

### OpenFOAM-10 Implementation  
- **Status**: ✅ Fully tested and validated
- **Features**: Complete CCM with ecVars system and full feature set
- **Testing**: ✅ Confirmed operational

### OpenFOAM-13 Implementation
- **Status**: ✅ Latest version support  
- **Features**: Complete CCM with ecVars system and OF-13 API compatibility
- **Testing**: On-going validation

## Error Control Variables (ecVars) Implementation

### Implementation Status
The ecVars system provides flexible variable encoding for CCM acceleration, allowing dynamic selection of important variables for thermochemical state grouping.

### System Design
- **Default Initialization**: ecVars initialized with non-species variables from principalVars (e.g., T, P, phi, chi)
- **EC Mode**: When enabled, ecVars evolves by adding/removing important species variables
- **Adaptive Evolution**: Variables added when importance > tolerance (automatically calculated as 1/nSlice)
- **Dynamic Management**: ecVars updated every `updateFreq` solver steps based on std analysis

### Key Implementation Details
1. **Initialization**: Constructor sets ecVars_ with non-species principalVars (line 314 in CCMchemistryModel.C)
2. **Tolerance Calculation**: Always `1.0/nSlice_` where `nSlice` is user-configurable  
3. **Variable Management**: Controlled by `numECVarsToAdd` and `numECVarsToRemove` parameters

### Configuration
```cpp
ecMode
{
    enabled          false;      // Enable adaptive variable selection
    numECVarsToAdd   3;         // Variables to add when expansion needed
    numECVarsToRemove 1;        // Variables to remove when optimization needed  
    updateFreq       5;         // Check every N solver steps
}
```


## Multi-Version Release Architecture

### Release Philosophy
This EC-CCM repository represents a mature, multi-version release supporting three major OpenFOAM versions (7, 10, 13) with consistent CCM functionality across all platforms.

### Version Compatibility Matrix
| OpenFOAM Version | Status | Features | Test Status |
|------------------|--------|----------|-------------|
| OF-7 | ✅ Stable | Full CCM with ecVars | On-going |
| OF-10 | ✅ Reference | Full CCM with ecVars + Complete Feature Set | ✅ Validated |  
| OF-13 | ✅ Latest | Full CCM with ecVars + Latest API | On-going |

### Shared Components
All versions implement identical core algorithms:
- `distributeReactionEntry()` - Cell grouping and MPI distribution
- `updateReactionRate()` - Chemistry computation with `getRRGivenYTP`
- `distributeReactionRate()` - Reaction rate distribution  
- `solve()` - Main chemistry solver integration

### Version-Specific Adaptations
- **OF-7**: Template syntax fixes (`this->` prefix, `template` keyword)
- **OF-10**: Reference implementation with full feature set
- **OF-13**: API compatibility updates for latest OpenFOAM

## Known Limitations by Version

### Testing Status by Version

### OpenFOAM-7 Status ✅
- Full CCM functionality including J/Phi/Chi calculations now operational
- Complete ecVars system implementation
- Testing: On-going validation to confirm all features

### OpenFOAM-10 Status ✅  
- All features fully tested and validated
- Complete ecVars system implementation
- Full J/Phi/Chi calculation support confirmed

### OpenFOAM-13 Status ✅
- Full CCM functionality implemented with OF-13 API compatibility
- Complete ecVars system implementation  
- Testing: On-going validation to confirm all features

## Recent Improvements

### Enhanced Debug Information (Latest)
**Implementation**: Updated species ranking output in `CCMDebug.C` for all versions
- **Change**: Modified debug output to clarify that only species with `norm_std > 1e-10` are displayed
- **Benefit**: Clearer understanding of displayed species selection criteria during debugging
- **Files Modified**: `OpenFOAM-7/code/helpers/CCMDebug.C:99`, `OpenFOAM-10/code/helpers/CCMDebug.C:99`, `OpenFOAM-13/code/helpers/CCMDebug.C:99`

### Optimized Flow Control (Latest)
**Implementation**: Improved procedure optimization after species removal in `CCMchemistryModel.C` for all versions  
- **Change**: Moved `distributeReactionEntry()` and `ecOkay` assignment inside the species removal condition block
- **Benefit**: Avoids unnecessary computation when no species are removed due to safety limits
- **Performance**: Reduces computational overhead in error control variable management
- **Files Modified**: All three versions' main chemistry model files at approximately line 1300-1310

## Quick Start Guide

### 1. Choose OpenFOAM Version
Select your target OpenFOAM version (7, 10, or 13) and activate the environment:
```bash
# Example for OpenFOAM-10
source /home/zhouyuchen/OpenFOAM/OpenFOAM-10/OpenFOAM-10/etc/bashrc
```

### 2. Compile CCM Library  
```bash
# Navigate to version-specific code directory
cd OpenFOAM-10/code/
wclean && wmake -j
```

### 3. Run Test Case
```bash
# Use version-specific test case directory
cd OpenFOAM-7/Sandia/    # For OF-7
cd OpenFOAM-10/Sandia/   # For OF-10  
cd OpenFOAM-13/Sandia/   # For OF-13

# Single core test
timeout 50s reactingFoam >& log.reactingFoam

# Or parallel test  
decomposePar -force
mpirun -n 6 reactingFoam -parallel >& log.reactingFoam&
```

## Running Procedure
1. Choose OpenFOAM version: Activate corresponding OpenFOAM environment (7, 10, or 13)
2. Compile code: Run `wmake -j` in respective version's `code/` directory  
3. Run test case: Enter version-specific test case directory (e.g., `OpenFOAM-10/Sandia/`) and run `reactingFoam`