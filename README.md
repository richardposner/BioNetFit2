# BioNetFit2

# New Features in the .conf file

The implementation for multiple models and model checking is complete for the GA, PSO, and DE. So they can be used as templates to update the SA algorithm.
I have also updated the model checking: now the user can provide the exact time point where the logical operation must be done, and the user can also choose which models must be compared.

Now constraints should be provided in .conf files in the following format:

      constraint=Param1[>,>=,==,<=,<]Param2 MODEL1 MODEL2 TIME1 TIME2

For example:

      constraint=RLbonds>=pR 0 1 0 60

where I'm comparing RLbonds from model 0 at time 0 versus pR from model 1 at time 60, and I'm checking if RLbonds is greater then or equal to pR



Another new option was included so the user can specify the constraint weight, which can be any value between 0 and 1
A very small number (i.e 0.1) = fit values matter the most
A very large number (i.e 0.9) = constraints matter the most
A balanced value (i.e 0.5) = both fit and contrsaint values are equally important

For example:

      constraint_weight=0.5

