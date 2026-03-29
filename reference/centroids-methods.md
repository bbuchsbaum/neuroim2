# Return a matrix of centroids of an object

Return a matrix of centroids of an object

## Usage

``` r
centroids(x, ...)

# S4 method for class 'ClusteredNeuroVec'
centroids(x, type = c("center_of_mass", "medoid"))

# S4 method for class 'ClusteredNeuroVol'
centroids(x, type = c("center_of_mass", "medoid"))
```

## Arguments

- x:

  an object with multiple centroids (e.g. a `ClusteredNeuroVol`)

- ...:

  extra args

- type:

  the type of center of mass: one of "center_of_mass" or "medoid"

## Value

A numeric `matrix` where each row represents the coordinates of a
centroid.

A matrix of coordinates where each row represents the centroid of a
cluster.

## Details

For \`type = "center_of_mass"\`, returns arithmetic mean coordinates;
for \`"medoid"\`, returns the most central point.
