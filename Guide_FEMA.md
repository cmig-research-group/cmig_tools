This section of the documentation is a guide on how to use FEMA for non-ABCD data. This is also meant to serve as a practical guide to explain how FEMA expects the data to be and what different things in FEMA mean.

The heart of FEMA is the ``FEMA_fit`` function that fits linear mixed effects (LME) model(s) on the input data. This is the function that you need to call for fitting LMEs.

# Notation:

``n``: observations  
``p``: predictors (fixed effects or *X* variables)  
``v``: voxels / vertices / ROIs / outcome variable (*y* variable)  
``c``: contrasts  
``r``: random effects  
``q``: number of unique observations ``r``: number of random effects  

# Inputs to ``FEMA_fit``:

``X``: numeric type ``n x p`` design matrix with fixed effects or explanatory variables  
``iid``: cell type which is ``n x 1`` indicating an ID for each observation; for example, same subjects will have the same ``iid``  
``eid``: cell type which is ``n x 1`` indicating an event ID for each observation; for example, all first visits across observations will have the same ``eid``; this is not currently used for analysis  
``fid``: cell type which is ``n x 1`` indicating the family ID for each observation; for example, siblings will (typically) have the same ``fid``  
``agevec``: numeric type which is ``n x 1`` indicating the age of each observation; this is not currently used for analysis  
``ymat``: numeric type ``n x v`` with outcome variables

# Additional inputs to ``FEMA_fit``:

``niter``: can be left as default = 1  
``contrasts``: vector of weights for testing hypothesis; can be left empty; see below for more details  
``nbins``: number that defines how many bins will be created; can be left as default = 20  
``pihatmat``: numeric ``q x q`` genetic relationship matrix; should be in the same order as ``unique(iid, 'stable')``; see below for more details  

Once these set of inputs are defined, ``FEMA_fit`` can be called as: ``FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, contrasts, nbins, pihatmat)`` with various outputs defined (see more below).  

# Optional additional inputs to ``FEMA_fit``:

``RandomEffects``:this controls which random effects will be accounted for (details in the next section)  
``nperms`` number of permutations; if greater than 1, permuted effects will also be output  
``FixedEstType``: default is ``GLS``; change to ``OLS`` if no random effects in the data  
``RandomEstType``: default is ``MoM``; change to ``ML`` for slower maximum likelihood-based estimation of variance components for the random effects  
``SingleOrDouble``: controls the numerical precision of computation; default is ``double``; can be set to ``single`` which can help with reducing memory requirement  
``PermType``: controls which type of permutation is performed; default is ``wildbootstrap``; can be alternativelt set to ``wildbootstrap-nn`` which is the non-null version  

# Random Effects:

To use FEMA, you need to define the random effects. Examples of random effects could be repeated measurements, family-like design (with or without repeated measurements), genetic relationship between subjects, etc. FEMA supports a large number of random effects:  

-   ``F``: family relatedness (related to the input ``FID``)
-   ``S``: subject - required for longitudinal analyses (related to the input ``IID``)
-   ``E``: error - always required
-   ``A``: additive genetic relatedness
-   ``D``: dominant genetic relatedness - square of A
-   ``M``: maternal effect - effect of having same mother
-   ``P``: paternal effect - effect of having same father
-   ``H``: home effect - effect of living at the same address
-   ``T``: twin effect - effect of having the same pregnancy ID

The most important random effect is the ``F`` or family effect: this defines the top level of clustering in the data. Therefore, observations are clustered within families. In case of simple longitudinal data (where there is no family structure), ``FID`` can be set to be the same as ``IID``.

The random effects are specified as a cell type vector:

-   If there is only family structure in the data, set ``RandomEffects = {'F', 'E'}`` (while setting the ``IID = FID``)
-   If there is only repeated measurement in the data, set ``RandomEffects = {'S', 'E'}`` (while setting the ``FID = IID``)
-   If there is family structure and repeated measurements, set ``RandomEffects = {'F', 'S', 'E'}``
-   If there is no family structure and no repeated measurements (i.e., your data satisfies the iid assumption), you don't really need mixed effects model and can just do linear regression. If you would like to fit these models in FEMA, you can set ``FID`` to unique values per observation, set ``IID = FID``, set ``RandomEffects = {'E'}``, and additionally pass, as input, ``FixedEstType`` which is equal to ``OLS``. This basically tells FEMA to fit a linear regression model instead of linear mixed effects model

If you have access to genetic data, you can create a genetic relatedness matrix (GRM), which can be calculated as the correlation coefficient of the standardized genotyping matrix. This GRM or ``pihatmat`` should be ordered as ``unique(iid, 'stable')`` and is ``q x q`` in size. The corresponding random effect can be specified as ``A``. For example, ``RandomEffects = {'F', 'A', 'E'}``, ``RandomEffects = {'F', 'A', 'S', 'E'}``.

Another thing to note is that family is just what we call the top level of clustering - for example, you could have subjects embedded within scanning sites (i.e., scanning site could be a random effect). In this case, ``FID`` can be set to scanner ID.

For other random effects beyond ``F``, ``A``, ``S``, and ``E``, additional optional inputs need to be passed. Specifically:  

- ``MotherID`` can be specified to model ``M`` or maternal effect
- ``FatherID`` can be specified to model ``F`` or paternal effect
- ``PregID`` can be specified to model ``T`` or twin effect
- ``HomeID`` can be specified to model ``H`` or home effect

One or more of these four IDs can be passed as inputs to ``FEMA_fit`` using the optional name-value pair combinations and the appropriate letter corresponding to random effect can be specified in the ``RandomEffects``.

Finally, note that dominant genetic relatedness or ``D`` can be modeled directly (i.e., no additional inputs beyond ``pihatmat`` is required; only the ``RandomEffects`` needs to include ``D``).

# Putting it all together

The call to ``FEMA_fit`` is ``FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, contrasts, nbins, pihatmat)`` along with any combination of optional additional inputs. The optional inputs are passed as name-value pairs. For example:

-   Calling ``FEMA_fit`` while specifying random effects: ``FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, contrasts, nbins, pihatmat, 'RandomEffects', {'F', 'A', 'S', 'E'})``
-   Calling ``FEMA_fit`` while specifying number of permutations to be 20: ``FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, contrasts, nbins, pihatmat, 'nperms', 20)``
-   Calling ``FEMA_fit`` while specifying numerical precision and random effects: ``FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, contrasts, nbins, pihatmat, 'RandomEffects', {'F', 'A', 'S', 'E'}, 'SingleOrDouble', 'single')``
-   Calling ``FEMA_fit`` while specifying pregnancy ID for modeling twin effect: ``FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, contrasts, nbins, pihatmat, 'RandomEffects', {'F', 'A', 'S', 'T', 'E'}, 'PregID', pregID)``, where ``pregID`` variable contains the list of pregnancy IDs for every observation in the data
-   Calling ``FEMA_fit`` while wanting to model the dominant genetic effect but not the additive genetic effect: ``FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, contrasts, nbins, pihatmat, 'RandomEffects', {'F', 'D', 'S', 'E'})``

# Explanation of main outputs:

``beta_hat``: numeric type ``c+p x v`` estimated beta coefficients; if contrasts are specified, the first ``c`` rows are the estimates for the contrasts  
``beta_se``: numeric type ``c+p x v`` estimated standard error for the beta coefficients; if contrasts are specified, the first ``c`` rows are the standard errors for the contrasts  
``zmat``: numeric type ``c+p x v`` estimated *Z* statistics (beta/SE); if contrasts are specified, the first ``c`` rows are the *Z* statistics for contrasts  
``logpmat``: numeric type ``c+p x v`` log10 *p* values; if contrasts are specified, the first ``c`` rows are the log10 *p* values for contrasts; note that the *p* values are signed, where the sign reflects the sign of the *Z* statistics  
``sig2tvec``: numeric type ``1 x v`` total residual variance for each ``v``  
``sig2mat``: numeric type ``r x v`` normalized estimates of variance components for each random effect in the order in which ``RandomEffects`` are defined  
The same variables with ``_perm`` in their names are estimates for the permuted data  

# Quick note on contrasts:

The contrast vector is numeric ``c x p`` vector/matrix where each row is a contrast that should be evaluated. A contrast can be used to define a hypothesis test. For example, if there are two *X* variables, then a contrast vector of ``[1 -1]`` tests the hypothesis that the difference between the two slopes (or beta coefficients of each of the two variables) is zero. A contrast vector ``[1 0]`` test the null hypothesis that the first coefficient is zero. By extension, this implies that if no contrasts are specified, all the beta coefficients returned as part of ``beta_hat`` are tested against zero. The contrast vector doesn't have to only contain zeros and ones. For example, a contrast vector ``[2 -1]`` would test the null hypothesis that the difference between two times beta coefficient 1 and beta coefficient is zero. A lot more can be said about the contrasts but the key to bear in mind is that contrasts are basically translated as ``cB`` or contrast vector multiplied by the beta coefficient vector. An important thing to bear in mind here is that in case of categorical variables, you should define your contrast carefully. If there are multiple dummy-coded (or multi-level) categorical *X* variables, the interpretation of the intercept is no longer the mean for reference level of one categorical variable; rather the intercept will then represent the mean across all reference levels for all categorical variables. This is, of course, not a FEMA-specific thing but has to do with the design matrix in general.

# Where to find example usage:

A useful starting point to get familiar with calling ``FEMA_fit`` is to look at the script ``FEMA_run_on_synthetic_data``; this script creates synthetic data and then performs the estimation by invoking ``FEMA_fit``. Another example is demonstrated in ``cmig_tools_utils/matlab/doFEMA_tests`` where the script creates simulated data for various combinations of random effects and then performs estimation. Additionally, in this script, for several combinations of random effects, we also compare the output using a standard LME solver. The simulations and analyses reported in the FEMA manuscript can be found [here](http://github.com/parekhpravesh/FEMA).

# How to get help:

Please feel free to open an issue on [GitHub](http://github.com/cmig-research-group/cmig_tools/issues) and we will be happy to help!

# Citation:

Parekh, P., Fan, C. C., Frei, O., Palmer, C. E., Smith, D. M., Makowski, C., Iversen, J. R., Pecheva, D., Holland, D., Loughnan, R., Nedelec, P., Thompson, W. K., Hagler, D. J. Jr, Andreassen, O. A., Jernigan, T. L., Nichols, T. E., & Dale, A. M. (2024). FEMA: Fast and efficient mixed-effects algorithm for large sample whole-brain imaging data. Human Brain Mapping, 45(2), e26579. <https://doi.org/10.1002/hbm.26579>
