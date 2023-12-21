Tree/Structure of the code:

top.m
?
?---start_pause.m
?
?---finish.m
?
?---display_rates.m
?
?---rates.m
??
??--- electron_part_balance.m
?
?---colour.m
?
?---rates_he.m, rates_argon.m, rates_he_n2.m, cr_model.m


RUNNING THE GLOBAL MODEL:
Copy the above files to your matlab local directory, open matlab and run top.m. 
You may need to change the "INPUT PARAMETERS SECTION" of top.m. Do not change the rest of the file!

GAS/CHEMISTRY:
The gas type can be selected in the "input parameters section" of top.m
To add a new gas to the global model, create a new rates_xxx.m file using 
rates_argon.m as an example. Do not forget to place the new file in the shared folder.
Modify top.m to accept a new gas option.   

CHANGES:
If you change top.m (apart from the input parameters section), finish.m, 
display_rates.m, rates.m, electron_part_balance.m