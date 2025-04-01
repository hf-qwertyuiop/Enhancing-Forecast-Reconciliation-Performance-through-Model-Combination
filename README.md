# Enhancing Forecast Reconciliation Performance through Model Combination
## Introduction
In result1.pdf, we report performance on different metrics, at different levels.

In result2.pdf, we take the labour data set as an example, select the best lambda with different metrics, and check its performance at each metric and level.

In result3.jpg, we show the win rate of our method on SMAPE (the main metric we report).

In E2E_accuracy_per_level.xlsx, we report the result of end-to-end approach at all hierarchical levels.

#### US macroeconomics series 1
The following code will make predictions from 20 training samples
```
python _main_make_predictions_for_first_macro_data.py
```
The output of forecasting accuracies in terms of APE and SIS at h=1,...,8 is 
