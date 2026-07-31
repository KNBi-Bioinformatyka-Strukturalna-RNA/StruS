# StruS


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

- Python 3.11+
- g++ compilator

### Installing

A step by step series of examples that tell you how to get a development env running.

#### [Linux]

Clone the repo and create virtual environment.

```bash
git clone https://github.com/KNBi-Bioinformatyka-Strukturalna-RNA/rna-model-error-detector.git na-model-error-detector
cd na-model-error-detector/StruS
python -m venv .venv
source .venv/bin/activate
(.venv) pip install -r requirements.txt
```

#### [Windows]

Clone the repo and create virtual environment.

```pwsh
git clone https://github.com/KNBi-Bioinformatyka-Strukturalna-RNA/rna-model-error-detector.git na-model-error-detector
cd na-model-error-detector/StruS
python -m venv .venv
.venv\Scripts\activate
(.venv) pip install -r requirements.txt
```

## Running the tool

### Single prediction RTBS:

```bash
(.venv) python StruS.py RTBS target.pdb prediction.pdb
```

### Multiple predictions structRMSD:

```bash
(.venv) python StruS.py structRMSD target.pdb -p predictions
```

### Run both tools with different output folder:

```bash
(.venv) python StruS.py target.pdb prediction.pdb -o results
```
