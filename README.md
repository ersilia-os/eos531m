# Human Liver Microsome Half-Life (PredHLM)

Predicts metabolic half-life in human liver microsomes from a SMILES, as log10 minutes and minutes, plus the spread of a five-model XGBoost ensemble on RDKit 2D descriptors and an applicability-domain score against the training set. Trained on 11,790 compounds from in-house Korea Chemical Bank assays and ChEMBL (test RMSE 0.507, R2 0.43). Phase I metabolism only, not prospectively validated; a decision-support tool, not a replacement for the assay.

This model was incorporated on 2026-10-15.


## Information
### Identifiers
- **Ersilia Identifier:** `eos531m`
- **Slug:** `hlm-half-life`

### Domain
- **Task:** `Annotation`
- **Subtask:** `Property calculation or prediction`
- **Biomedical Area:** `ADMET`
- **Target Organism:** `Homo sapiens`
- **Tags:** `Half-life`, `Microsomal stability`, `ADME`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `4`
- **Output Consistency:** `Fixed`
- **Interpretation:** Predicted human liver microsome half-life as log10 minutes and minutes, with ensemble spread and training-set similarity score.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| pred_log10_half_life | float | high | Predicted log10 metabolic half-life (minutes) in human liver microsomes (ensemble mean of 5 XGBoost models) |
| pred_std | float |  | Standard deviation of the 5 fold-model predictions in log10 units (uncertainty) |
| pred_half_life_min | float | high | Predicted metabolic half-life in human liver microsomes in minutes (10 to the power of the log10 prediction) |
| sdc | float | high | Applicability-domain score (sum of distance-weighted contributions from ECFP4-similar training compounds); higher means better supported by the training set |


### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos531m.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos531m.zip)

### Resource Consumption
- **Model Size (Mb):** `34`
- **Environment Size (Mb):** `1506`


### References
- **Source Code**: [https://github.com/joshua-416/PredHLM](https://github.com/joshua-416/PredHLM)
- **Publication**: [https://doi.org/10.64898/2026.07.02.736062](https://doi.org/10.64898/2026.07.02.736062)
- **Publication Type:** `Preprint`
- **Publication Year:** `2026`
- **Ersilia Contributor:** [TiagoJanela](https://github.com/TiagoJanela)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [MIT](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos531m
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos531m
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
