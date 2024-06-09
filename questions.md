> One important issue is with electric and magnetic polarization modes, from what I see you order first all electric and then all magnetic, but in the standard now those are alternating. Honestly, I am not fully aware of the advantages, it came from the definition in treams, do not know if there are advantages for you in your way of storing.

- No problem
- Does not seem very intuitive

> As I see, in your parameter sweep you plan to save in the same file only wavelength variation, and you do not have any method parameters that change for different wavelengths, so this is all in order.  

I think this was decided at the workshop

> But the sweep is somehow inner, not the outer dimension, this can cause problems. 

To clarify

> In computation, now the approach is to have computation/files/ group that includes the source files needed to repeat the computation, and computation/method_parameters group for all the additional parameters that are used by the software.  


'method','EBCM',...
    'software','SMARTIES',...
    'version','1.1',...
    'unit','nm', ...
    'Lmax', globalN, ...
    'Ntheta', globalnNbTheta, ...
    'accuracy', accuracy
    
Should add "script", currently in comments as

computation/files 

and group the others under

method_parameters

> The name of the material is easier to search in the database if specified as an attribute, so the group itself can be named however you like, but the attribute "name" which has the value "Ag" or "Silver"  or better both, would be the place to search for the material name. The reference to permittivity values(as a link to refractiveindex.info etc.) can be in the attribute "reference" of the group, could add also the interpolation method name in the attribute "interpolation" (those are some new aspects that came after discussion in the workshop).  

Add "name" attribute to group as keywords Ag, Silver
Add "reference" attribute with DOI
Add "interpolation" attribute


> Similarly in geometry, the group name should match the group name in materials/, and have the attribute  "shape" to set the name of the shape as you have it, datasets like radius etc are part of the group. 

- group name should match the materials group
- radius etc. are part of the group

> You also specify shape as dataset with list of letters, this we did not have before, I do not know whether this brings extra advantages.  

Don't know what this means

> You are specifying in the name of the shape  that the spheroid is oblate or prolate, which we did not do before, but this is searchable under spheroid so is fine. 

OK

> "unit" can be added as an attribute to this group, if you do not provide a mesh. 

add "unit" attribute to group


