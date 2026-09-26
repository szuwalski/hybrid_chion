# Methods revision: style and scope notes

## Style adopted from the reference

The Methods now open with the purpose and overall operation of the joint
model. After the input data, the population section introduces the state
variables and annual cycle, then follows survival and fishing, growth,
recruitment, terminal molt, and survival to the next census. The observation
model follows the biological transitions; the objective, regularization, and
estimation follow the observation model. This uses the reference's connected,
process-first explanation and introduces equations after explaining what
they do. Numerical transforms and most prior constants have moved to
Supplementary methods within manuscript.Rmd rather than being discarded.

## Where the reference's approach needs adaptation

- **Separate-stock framing:** the reference introduces largely separate male
  population models. This manuscript needs an early explanation of six
  components, females, hybrid constraints, and shared classification error.
  Presenting these as independent stock fits would obscure the central model.
- **One annual sequence for every population:** a common narrative is useful,
  but the current implementation uses a split mortality interval and fishing
  for parental males, and full annual mortality before growth for females and
  hybrid males. The revised equations explicitly preserve that difference.
- **Equation-by-equation detail:** the reference often interrupts the process
  narrative with several scalar equations and parameter lists. The revision
  uses shared vector transitions and defers numerical regularization details
  to the supplement. It keeps definitions near the equations.
- **Terminology and notation:** the reference alternates size/time symbols and
  temporary abundance notation. The revision consistently uses g, t, z, and
  I/M, and distinguishes captured, retained, and surviving crab. It also
  distinguishes biological hybrids from recorded labels. Classification
  probability is not itself a hybridization rate or mating probability.
- **Borrowed technical claims:** the reference's externally estimated growth,
  time-invariant fishery ogives, king-crab molting, second-difference mortality
  penalties, and ADMB implementation do not describe this model. Here growth
  is jointly estimated, fishery ogives vary annually, mortality deviations
  have independent penalties, and estimation uses RTMB. Only the explanatory
  structure was borrowed. The reference's displayed likelihood and growth
  formulas were not copied; equations follow the current implementation.
- **Strength of claims:** phrases such as “best available information” need
  support, and claims of heatwave causation or directly measured hybridization
  do not follow from population accounting alone. The opening describes what
  the model estimates without asserting a mechanism it has established.
- **Diagnostics as completed actions:** the reference reports 100 jittered
  starts and retrospective analyses. Those cannot be adopted as prose unless
  performed for this model. The current fit's small gradient and positive
  Hessian do not erase its false-convergence flag; final diagnostics should
  explicitly state the checks actually completed and their limitations.

## Substantive issues to resolve before submission

1. **Hybrid recruitment scope:** the prior draft deliberately omitted the
   internally fitted parental-abundance effect on hybrid recruitment. That
   choice and its provenance warning are preserved. The current JMR figures
   nevertheless come from a fit containing that effect. A publication cannot
   leave the hybrid entry flux incompletely specified while presenting those
   figures as results of the stated model. Describe the fitted relationship,
   or refit the intended manuscript model and replace the figures. This is
   an analysis/specification decision, not a stylistic edit.
2. **New snow-male bin:** the Methods now state the shared 27.5 mm minimum and
   describe how the added snow-male observations were handled. Retaining the
   original survey CVs, extrapolating selectivity, assuming zero terminal
   molt in years without the smallest-class maturity input, and assigning
   zero catch composition at the new bin are assumptions requiring disclosure
   and sensitivity assessment. They should not disappear for smoother prose.
3. **Incomplete source and analysis sections:** survey details and citations,
   component coverage, composition-weight derivation, environmental analyses,
   and hybridization analyses still need substantive completion. Placeholders
   remain rather than importing the reference's species, covariates, GAM
   settings, management projections, or diagnostics without verification.
4. **Results and editorial material:** the scope comments, bracketed requests,
   and figure provenance note are useful while drafting, but should be
   resolved before submission. Exact fitted gradients and Hessian results
   belong in the final diagnostics/results or supplement; the Methods should
   describe how numerical adequacy was assessed.

The reference manuscript, model code, fitted objects, and non-Methods prose
were not revised for this editorial task. A supplementary methods subsection
was added to retain technical details moved out of the main Methods.
