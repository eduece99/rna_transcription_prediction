# AGENTS.md

## Project overview
This workspace contains a synthetic computational biology exercise for selecting durable liver-expression UTR designs. The goal is to predict which untested constructs are most promising and to propose new candidate designs based on the provided data.

Primary materials:
- `brief.md` — task summary and required outputs
- `primer.md` — biological background and interpretation guidance
- `data_dictionary.md` — schema definitions for all CSV files
- `utr_library.csv` — all 3,150 construct sequences and metadata
- `invitro_screen.csv` — cell-line screen with reporter and RNA measurements
- `invivo_liver.csv` — mouse liver RNA data for 60 tested constructs
- `prediction_targets.csv` — target constructs for prediction

## Working rules for agents
- Keep all work inside this repository and avoid creating machine-specific paths.
- Treat the supplied CSV files as read-only inputs unless explicitly asked to transform them in a reproducible script.
- Prefer transparent, interpretable modeling approaches over opaque black-box pipelines when the data and task are small enough.
- Document assumptions clearly, especially around synthesis-batch effects, animal/cage variation, and assay noise.
- Preserve the synthetic nature of the data and do not claim biological conclusions beyond what the evidence supports.
- Do not share or publish the data files outside the workspace.

## Expected deliverables
The final solution should produce these artifacts in the repo root or an organized subfolder:
1. `predictions.csv` with columns:
   `construct_id,predicted_day28,lower80,upper80`
2. `designs.csv` with columns:
   `design_id,utr5,utr3,predicted_day28,lower80,upper80,rationale`
3. `report.md` written for bench scientists and no more than two pages
4. `code/` containing reproducible analysis scripts, environment/dependency metadata, and fixed seeds for any randomized steps

## Recommended workflow
1. Read `brief.md`, `primer.md`, and `data_dictionary.md` before modeling.
2. Inspect the structure and relationships in the CSV files to identify key signals:
   - UTR sequence features relevant to translation and degradation
   - construct family/parent/variant relationships
   - batch structure and replicate patterns
   - observed variation across in vitro vs in vivo data
3. Build a small, explainable predictive workflow that links sequence features to day-28 liver RNA.
4. Validate with transparent baselines, such as a simple regression or a naive benchmark, and compare against the final approach.
5. Produce uncertainty intervals that reflect the actual uncertainty of a new five-mouse study mean rather than an overly optimistic confidence band.
6. Keep all code reproducible and executable from a single documented command.

## Reproducibility expectations
- Use fixed random seeds anywhere randomness is involved.
- Include a single command that rebuilds all required outputs.
- Record runtime and language version in the project code or generated docs.
- Keep code dependencies explicit via a requirements file, environment file, or equivalent.
- Use relative paths only.

## Reporting standards
- Lead with a recommendation, not methodology alone.
- Explain what was tried and why the final approach was chosen over alternatives.
- State how confident the predictions are and what evidence supports that confidence.
- Describe what would change the recommendation.
- Name the next measurements that would answer the key remaining questions.

## High-priority interpretation notes
- The synthetic benchmark is designed to reward reasoning and honest uncertainty over perfect accuracy.
- The 5'UTR mainly affects translation efficiency; the 3'UTR mainly affects transcript abundance and decay.
- miR-122 seed motif `ACACUCC` and AU-rich motifs such as `AUUUA` are strong biological features to consider.
- Mouse data is measured at day 3 and day 28 in separate cohorts, which should be accounted for when estimating uncertainty and consistency.
- 80% predictive intervals should reflect uncertainty in the mean of a future five-mouse study using newly prepared material.

## Suggested repository conventions
- Keep scripts under `code/` with a clear entry point, such as `code/run_analysis.py`.
- Use concise docstrings and comments explaining non-obvious modeling choices.
- Preserve one source of truth for metadata and output generation.
- Avoid large, opaque notebooks unless they are clearly reproducible and well documented.

## Completion checklist
Before finishing, confirm that the workspace contains:
- a clear analysis pipeline
- reproducible outputs
- honest uncertainty estimates
- a benchmark comparison
- a concise scientist-facing report
- no hardcoded local paths or unexplained assumptions
