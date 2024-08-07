function [colnames_model X] = FEMA_DEAP_job_data(fstem_imaging,fname_design,dirname_tabulated,dirname_imaging)

  if nargin < 4
    logging('Usage: FEMA_DEAP_wrapper(fstem_imaging,fname_design,dirname_tabulated,dirname_imaging)');
    error('Incorrect number of input arguments')
  end

  % Ben: is this needed?
  rng shuffle % Set random number generator so different every time -- should allow for option to control this

  % Create analysis job
  [colnames_model X] = FEMA_DEAP_wrapper_submit(fstem_imaging,fname_design,dirname_tabulated,dirname_imaging,'voxel');

  return