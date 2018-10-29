[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 50
[]

[Problem]
  solve = false
[]

[AuxVariables]
  [u_aux]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[ICs]
  [u_aux]
    type = VolumeWeightedWeibullIC
    legacy_generator = false
    variable = u_aux
    reference_volume = 100
    weibull_modulus = 15.0
    median = 1.0
  []
[]

[VectorPostprocessors]
  [./histo]
    type = VolumeHistogram
    variable = u_aux
    min_value = 0
    max_value = 4
    bin_number = 80
    execute_on = initial
    outputs = initial
  [../]
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  [./initial]
    type = CSV
    execute_on = initial
  [../]
[]
