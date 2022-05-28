# BLAMlab

This is a repository for a motor learning study in BLAM Lab.

## Instructions for data path in analyzation
Data was stored under the **"subject"** folder, sorted by different subject. Under the folder of each subject, data was sorted with different tasks. 
| Abbr. | Description |
| ----------- | ----------- |
| avg | Average mapping |
| bi | Bimanual mapping |
| C | Spatial Working Memory (Corsi) |
| R | Mental Rotation |
| n | Not involving |
| d1 | day 1 |

**List of all the tasks**: 
avg-nC-nR, avg-C-nR, avg-nC-R, avg-C-R, bi-nC-nR, bi-C-nR, bi-nC-R, bi-C-R, d1_prac_avg1, d1_prac_bi1, d1_prac_bi2, d1_prac_bi3, d1_prac_bi4, d1_prac_bi5, d1_prac_bi6, d2_prac_bi7

- Example 1: "Data\subject\subject\subj05\avg-C-nR"
- Example 2: "Data\subject\subject\subj05\d1_prac_bi1"

## Description for the code and data
| File name | Description |
| ----------- | ----------- |
| taskOrder_02232022.m | Generating the order of 2x2 tasks for different subjects |
| taskOrder_02232022.txt | The order of 2x2 tasks for the 16 subjects |
| corsiOrder_02222022.m | Generating the tFile of differnt trial in each task |
| CorsiRotation_05112022.m | Main analyzation for the reaching direction error, reaching path length, movement duration |
| ThesisData.mat | Data file for running CorsiRotation_05112022.m |
| repeated_anova_0426.R | Repeated anova R file |
| rdata.csv | Data file for the Repeated anova in R |
| LearningCurve_05122022.m | Learning curve of bimanual mapping represented by normalized path length by one subject|
| LearningCurve.mat | Data file for running LearningCurve_05122022.m |
