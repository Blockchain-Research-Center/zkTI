# Crowdsourcing TI Methods

These codes and datasets are implementation of truth inference (TI) methods mentioned in the survey: [Truth inference in crowdsourcing: is the problem solved?](https://doi.org/10.14778/3055540.3055547). They are Python codes, and the datasets applied are in the `../datasets` folder.

Parts of the codes are borrowed from the open sourcing library: [crowd_truth_infer](https://github.com/zhydhkcws/crowd_truth_infer)

## Annotation

These annotations followes the base work. See the open source algorithm library above for more details.

- w: workers
- l: labels
- e: examples, or tasks
- lpd: example's predicted labels with its probability

## Usage

``` python
python3 $algorithm_script$ $answer_file$ $truth_file$
```

For example, 

``` python
python3 MV.py ../datasets/d_duck_identification/answer.csv ../datasets/d_duck_identification/truth.csv
```