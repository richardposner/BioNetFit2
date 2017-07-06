# BioNetFit2
The implementation for multiple models and model checking is complete for the GA, PSO, and DE. So they can be used as templates to update the SA algorithm.
I have also updated the model checking: now the user can provide the exact time point where the logical operation must be done, and the user can also choose which models must be compared.

Now constraints should be provided in .conf files in the following format:
      constraint=Param1[>,>=,==,<=,<]Param2 MODEL1 MODEL2 TIME1 TIME2

For example:
      constraint=RLbonds>=pR 0 1 0 60
where I'm comparint RLbonds from model 0 at time 0 versus pR from model 1 at time 60, and I'm checking if RLbonds is greater then or equal to pR
