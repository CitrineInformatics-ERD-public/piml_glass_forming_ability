# Physics-informed machine learning for glass forming ability

Supporting datasets and scripts for the paper "Evaluation of 
GlassNet for physics-informed machine learning of glass stability 
and glass forming ability", Allec et al. A preprint of the 
manuscript is available on arXiv.

Jupyter notebooks for reproducing the figures of the paper are in
the [`notebooks`](notebooks) folder. Supporting datasets are in the
[`data`](data) folder.

## Python environment setup

After activating your python environment that you will be using
to run these notebooks (e.g., `conda activate piml_gfa`), use pip
to install the requirements:

```
pip install -r requirements.txt
```

Note that installation of [`glasspy`](https://github.com/drcassar/glasspy) requires [`pytorch`](https://pytorch.org/) (for [`pytorch`](https://pytorch.org/)
installation, see [here](https://pytorch.org/get-started/locally/).

Also, different notebooks require different [`glasspy`](https://github.com/drcassar/glasspy) versions in order to
reproduce the results in the paper exactly. Please see the top cell of each
notebook to see if a specific version is recommended.

## Running the notebooks
All of the notebooks can be ran independently of each other; however,
certain lines are time consuming and can be reduced by saving data 
locally. Two examples of this are i) the loading of the GlassNet training
and test datasets and ii) the regression of the MYEGA equation for these 
datasets.

### 1. Loading GlassNet data
Currently, in each notebook, we use [`glasspy`](https://github.com/drcassar/glasspy) to load the GlassNet data
like this:

```
# Import GlassNet
from glasspy.predict.models import GlassNet

# Load the model
glassnet = GlassNet()

# Then load data
glassnet_test_df = glassnet.get_test_dataset()
glassnet_train_df = glassnet.get_training_dataset()
```

Instead, you can do this once on your machine, and then save these dataframes
as csv files and read them in later:

```
# Save dataframes to csv files
glassnet_test_df.to_csv('glassnet_test_df.csv')
glassnet_train_df.to_csv('glassnet_train_df.csv')

# Load dataframes from csv files
glassnet_test_df = pd.read_csv('glassnet_test_df.csv', index_col=0, header=[0, 1])
glassnet_train_df = pd.read_csv('glassnet_train_df.csv', index_col=0, header=[0, 1])
```

### 2. MYEGA regression
Currently, in any notebook that requires $\eta(T_l)$, we use GlassNet to regress
the MYEGA equation:

```
# Import GlassNet
from glasspy.predict.models import GlassNet

# Load the model
glassnet = GlassNet()

# Read training data from csv file (this can only be done after completing 1. above)
glassnet_train_df = pd.read_csv('glassnet_train_df.csv', index_col=0, header=[0, 1])

X = glassnet_train_df.elements
Tl = glassnet_train_df.property.Tl
visc_at_Tl = glassnet.predict_log10_viscosity( T = Tl, composition = X )
```

Similar to the example above, you can do this once and then save the viscosity in a
csv file to load later:

```
# Save series to csv file
visc_at_Tl.to_csv('viscosity_glassnet_train.csv')

# Load series from csv file
visc_at_Tl = pd.read_csv('viscosity_glassnet_train.csv', index_col=0)
```

## Other notes
- For the ternary plots, the only systems for which the true glass forming region 
can be plotted are Na2O-Fe2O3-P2O5 and Na2O-B2O3-SiO2. This is hard-coded. If you'd
like to add your own system, you will have to edit `plot_ternary` of `notebooks/utils.py`. 


