Phylogeny, Biogeography, and Ancestral State Reconstructon of Trichosternus series
================
Sean Perez
April 3, 2018

Overview
--------

I am constructing a phylogeny for a diverse clade of ground beetles (Coleoptera:Carabidae) known as the Trichosternus series.

![Alt text](/Users/seanperez/Desktop/Trich_Series/ESPM288_Final_Project_Phylogeny/Trichosternus%20vigorsi.png)

*Trichosternus vigorsi* copyright copyright Udo Schmidt.

Data Source
-----------

I will be using DNA sequence data from four genes obtained from specimens in my lab.

In order to infer the historical biogeography, I will download, clean, and aggreggate data from the Essig Museum Database which contains locality information for my species of interest.

Goals
-----

1.  Construct a phylogeny through parsimony, maximum likelihood, and bayesian inference.

2.  Reconstruct ancestral character of brooding, a behavior not known in other ground beetles (Carabidae).

3.  Use locality information from the Essig Museum database to map species ranges and answer hypotheses about historical biogeograph.

Downloading locality data for species of interest from gbif.
============================================================

Install the 'rgbif' package.

    ## ── Attaching packages ────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 2.2.1.9000     ✔ purrr   0.2.4     
    ## ✔ tibble  1.4.2          ✔ dplyr   0.7.4     
    ## ✔ tidyr   0.8.0          ✔ stringr 1.3.0     
    ## ✔ readr   1.1.1          ✔ forcats 0.3.0

    ## ── Conflicts ───────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ✖ dplyr::vars()   masks ggplot2::vars()

Find the relevant taxa.

``` r
(key <- name_suggest(q='Nurus', rank = 'genus')$key[1])
```

    ## [1] 4751178

A general overview of the data quantity contained in Gbif.

``` r
res <- occ_search(taxonKey = key, limit=100)
res
```

    ## Records found [839] 
    ## Records returned [100] 
    ## No. unique hierarchies [8] 
    ## No. media records [1] 
    ## No. facets [0] 
    ## Args [limit=100, offset=0, taxonKey=4751178, fields=all] 
    ## # A tibble: 100 x 79
    ##    name       key decimalLatitude decimalLongitude issues   datasetKey    
    ##    <chr>    <int>           <dbl>            <dbl> <chr>    <chr>         
    ##  1 Nurus   1.63e9           -26.8             153. cdround… 84a649ce-ff81…
    ##  2 Nurus … 1.45e9           -28.9             153. cdrep    0645ccdb-e001…
    ##  3 Nurus … 1.45e9           -28.8             153. cdrep    0645ccdb-e001…
    ##  4 Nurus … 1.45e9           -28.9             153. cdrep    0645ccdb-e001…
    ##  5 Nurus   1.06e9           -28.6             153. gass84,… 5d283bb6-64dd…
    ##  6 Nurus   1.06e9           -28.2             153. gass84,… 5d283bb6-64dd…
    ##  7 Nurus … 1.63e9           -28.8             153. gass84   a79c2b50-6c8a…
    ##  8 Nurus … 1.06e9           -21.2             149. gass84,… 5d283bb6-64dd…
    ##  9 Nurus   1.06e9           -21.1             149. gass84,… 5d283bb6-64dd…
    ## 10 Nurus … 1.06e9           -21.1             149. gass84,… 5d283bb6-64dd…
    ## # ... with 90 more rows, and 73 more variables: publishingOrgKey <chr>,
    ## #   publishingCountry <chr>, protocol <chr>, lastCrawled <chr>,
    ## #   lastParsed <chr>, crawlId <int>, extensions <chr>,
    ## #   basisOfRecord <chr>, individualCount <int>, taxonKey <int>,
    ## #   kingdomKey <int>, phylumKey <int>, classKey <int>, orderKey <int>,
    ## #   familyKey <int>, genusKey <int>, scientificName <chr>, kingdom <chr>,
    ## #   phylum <chr>, order <chr>, family <chr>, genus <chr>,
    ## #   genericName <chr>, taxonRank <chr>, year <int>, month <int>,
    ## #   day <int>, eventDate <chr>, modified <chr>, lastInterpreted <chr>,
    ## #   license <chr>, identifiers <chr>, facts <chr>, relations <chr>,
    ## #   geodeticDatum <chr>, class <chr>, countryCode <chr>, country <chr>,
    ## #   eventID <chr>, identifier <chr>, occurrenceStatus <chr>,
    ## #   recordedBy <chr>, institutionID <chr>, datasetName <chr>,
    ## #   eventTime <chr>, datasetID <chr>, gbifID <chr>, occurrenceID <chr>,
    ## #   speciesKey <int>, species <chr>, specificEpithet <chr>,
    ## #   coordinatePrecision <dbl>, coordinateUncertaintyInMeters <dbl>,
    ## #   stateProvince <chr>, http...unknown.org.classs <chr>, county <chr>,
    ## #   collectionCode <chr>, catalogNumber <chr>, vernacularName <chr>,
    ## #   institutionCode <chr>, ownerInstitutionCode <chr>,
    ## #   occurrenceRemarks <chr>, lifeStage <chr>, elevation <dbl>,
    ## #   references <chr>, preparations <chr>, rights <chr>, locality <chr>,
    ## #   accessRights <chr>, locationRemarks <chr>, verbatimLocality <chr>,
    ## #   identifiedBy <chr>, elevationAccuracy <dbl>

There are 839 records of the genus Nurus in this dataset, a reasonable amount to work with if all are imported. Setting limit to 1000 since the default is only 500.

``` r
Nurus_gbif <- occ_search(taxonKey = key, limit = 1000)
```

Checking for all unique species names to see if any erroneous taxa were imported.

``` r
uniqueNurus <- unique(Nurus_gbif$data$species)
uniqueNurus
```

    ##  [1] NA                 "Nurus atlas"      "Nurus brevis"    
    ##  [4] "Nurus medius"     "Nurus grandis"    "Nurus curtus"    
    ##  [7] "Nurus nox"        "Nurus niger"      "Nurus latipennis"
    ## [10] "Nurus imperialis" "Nurus rex"        "Nurus fortis"

Filtering out for just species and their respective localities.

``` r
Nurus_localities <- select(Nurus_gbif$data, "species", "decimalLatitude", "decimalLongitude", "key")

head(Nurus_localities)
```

    ## # A tibble: 6 x 4
    ##   species      decimalLatitude decimalLongitude        key
    ##   <chr>                  <dbl>            <dbl>      <int>
    ## 1 <NA>                   -26.8             153. 1632976095
    ## 2 Nurus atlas            -28.9             153. 1452756101
    ## 3 Nurus brevis           -28.8             153. 1452756089
    ## 4 Nurus atlas            -28.9             153. 1452756095
    ## 5 <NA>                   -28.6             153. 1060924928
    ## 6 <NA>                   -28.2             153. 1060922902

Counting the number of NAs in the species column.

``` r
nrow(Nurus_localities)
```

    ## [1] 839

``` r
sum(is.na(Nurus_localities$species))
```

    ## [1] 143

143 out of 839 records are 'NA' for species, potentially representing species not identified to species. These will be removed from the analysis.

``` r
Nurus_localities_omitted <- na.omit(Nurus_localities, cols = "species")
```

Checking for any leftover NAs and new number of rows.

``` r
sum(is.na(Nurus_localities_omitted$species))
```

    ## [1] 0

``` r
nrow(Nurus_localities_omitted)
```

    ## [1] 672

Adding records to NCBI
----------------------

No nucleotide information is currently available for species in the Trichosternus-complex. This project aims to broadly sample clades and make them available for public use.
