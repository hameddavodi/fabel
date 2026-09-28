# Machine learning: scikit-learn pipelines and PyTorch layers

fdatools' estimators follow the scikit-learn API (`fit`, `transform`,
`predict`), so they drop into a `Pipeline` and a `GridSearchCV` with no
wrappers. The optional `fdatools.nn` module adds PyTorch layers, so a smoothing
step or a basis evaluation can sit inside a neural network and be trained
with the rest of it.

Our task: tell boys from girls using only their height curves from the
Berkeley growth study (it ships with fdatools). Boys grow for longer and have a
later, stronger growth spurt, so the curves carry the answer.

## 1. The data as a table

scikit-learn wants one row per sample. Here a sample is a child, and the
columns are the heights at the 31 ages.

```python
import matplotlib.pyplot as plt
import numpy as np
import fdatools as fdt

growth = fdt.datasets.load_growth()
age = growth.age
X = np.hstack([growth.hgtm, growth.hgtf]).T            # (93, 31): 39 boys, 54 girls
y = np.r_[np.zeros(39), np.ones(54)]                   # 0 = boy, 1 = girl
X.shape, y.mean().round(2)
```

## 2. A pipeline: smooth, reduce, classify

The pipeline has three steps:

1. `Smoother` turns each row of raw heights into basis coefficients
   (penalised smoothing, as in `fdt.smooth`).
2. `FPCA` turns the coefficients into a few principal component scores.
3. `LogisticRegression` classifies the children from those scores.

Both fdatools steps must use the same basis, so we pass it to both. `FPCA` then
measures distances between curves in the right metric (the integral of the
squared difference), not just between coefficient vectors.

```python
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.pipeline import Pipeline

from fdatools.decomposition import FPCA
from fdatools.smoothing import Smoother

basis = fdt.BSpline(domain=(1.0, 18.0), n_basis=20, order=6)
pipe = Pipeline([
    ("smooth", Smoother(basis, t=age, lam=1e-2, penalty=3)),
    ("fpca", FPCA(n=3, basis=basis)),
    ("clf", LogisticRegression(max_iter=1000)),
])

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=0, stratify=y
)
pipe.fit(X_train, y_train)
round(pipe.score(X_test, y_test), 2)                   # accuracy on unseen children
```

## 3. Tuning with `GridSearchCV`

Every setting of every step is a parameter named `<step>__<setting>`, so
`GridSearchCV` can tune the smoothing level and the number of harmonics
together, with 5-fold cross-validation:

```python
search = GridSearchCV(
    pipe,
    {"fpca__n": [2, 3, 5], "smooth__lam": [1e-2, 1.0]},
    cv=5,
)
search.fit(X_train, y_train)

search.best_params_
round(search.best_score_, 2), round(search.score(X_test, y_test), 2)
```

The fitted steps stay available for inspection. For example, the share of
variance of each harmonic in the best model:

```python
best = search.best_estimator_
best.named_steps["fpca"].varprop.round(3)
```

`Smoother`, `FPCA`, `FCCA`, `Registrator`, `FRegress` and `PDA` all follow
this API.

## 4. PyTorch: a dataset of curves

The rest of this page needs the `torch` extra:
`uv add "fdatools[torch]"`.
`import fdatools` never imports PyTorch; `fdatools.nn` loads it on first use.

`FDataDataset` serves the curves of an `FData` to a PyTorch `DataLoader`.
Each item is one curve's coefficients (or its values on a grid, with `t=`)
together with its label.

```python
# requires: torch
import torch
from fdatools.nn import BasisLayer, FDataDataset, SmoothingLayer

torch.manual_seed(0)
curves = fdt.smooth(X.T, age, basis=basis, lam=1e-2, penalty=3).fd   # 93 curves
dataset = FDataDataset(curves, y)
loader = torch.utils.data.DataLoader(dataset, batch_size=16, shuffle=True)

coefs, label = dataset[0]
coefs.shape, len(dataset)                              # (20,) coefficients, 93 curves
```

## 5. A network that reads growth speed

`BasisLayer` evaluates a basis expansion: it maps coefficients to curve
values on a grid. With `deriv=1` it gives the growth speed. A linear layer on
top of the speed curve makes a small classifier, trained with an ordinary
PyTorch loop:

```python
# requires: torch
grid = np.linspace(1.0, 18.0, 69)
model = torch.nn.Sequential(
    BasisLayer(basis, grid, deriv=1),                  # coefficients -> speed at 69 ages
    torch.nn.Linear(69, 1, dtype=torch.float64),
)
optimizer = torch.optim.Adam(model.parameters(), lr=1e-2)
loss_fn = torch.nn.BCEWithLogitsLoss()

for epoch in range(60):
    for coefs, label in loader:
        optimizer.zero_grad()
        loss = loss_fn(model(coefs)[:, 0], label.to(torch.float64))
        loss.backward()
        optimizer.step()

with torch.no_grad():
    predicted = (model(dataset.data)[:, 0] > 0).numpy()
round(float((predicted == y).mean()), 2)               # training accuracy
```

The weights of the linear layer are a function of age: they show which part
of the speed curve the network uses.

```python
# requires: torch
weights = model[1].weight.detach().numpy()[0]
fig, ax = plt.subplots()
ax.plot(grid, weights)
ax.axhline(0.0, color="grey", linewidth=0.5)
ax.set(xlabel="age (years)", ylabel="weight", title="What the network looks at")
```

## 6. Learning the smoothing level

`SmoothingLayer` is penalised smoothing as a layer: it maps raw observations
to coefficients. With `trainable_lam=True` the smoothing parameter λ becomes
a model parameter (stored as log λ), and it is trained by gradient descent
together with everything else. The network now starts from the raw heights.

```python
# requires: torch
smoother = SmoothingLayer(basis, age, lam=1.0, penalty=3, trainable_lam=True)
network = torch.nn.Sequential(
    smoother,                                          # raw heights -> coefficients
    BasisLayer(basis, grid, deriv=1),                  # coefficients -> speed
    torch.nn.Linear(69, 1, dtype=torch.float64),
)
raw = torch.as_tensor(X)                               # (93, 31)
target = torch.as_tensor(y)
optimizer = torch.optim.Adam(network.parameters(), lr=1e-2)

for step in range(200):
    optimizer.zero_grad()
    loss = loss_fn(network(raw)[:, 0], target)
    loss.backward()
    optimizer.step()

round(smoother.lam, 2)                                 # lambda has moved away from 1.0
```

Gradients flow through the smoothing solve, the basis evaluation and the
classifier. On a GPU, move the model with `network.to("cuda")`; on Apple
silicon use `dtype=torch.float32` in the layers, since MPS has no float64.

## R equivalent

R's `fda` has no machine-learning integration. The closest is to build the
steps by hand:

```r
library(fda)
hgt      <- cbind(growth$hgtm, growth$hgtf)
sex      <- c(rep(0, 39), rep(1, 54))
basis    <- create.bspline.basis(c(1, 18), 20, norder = 6)
hgtfd    <- smooth.basis(growth$age, hgt, fdPar(basis, 3, 1e-2))$fd
scores   <- pca.fd(hgtfd, nharm = 3)$scores
fit      <- glm(sex ~ scores, family = binomial)
```

Tuning then needs your own cross-validation loop, and nothing is
differentiable end to end. In fdatools the same steps are one `Pipeline`, and
`fdatools.nn` makes them trainable layers.
