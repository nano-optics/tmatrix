
> One important issue is with electric and magnetic polarization modes, from what I see you order first all electric and then all magnetic, but in the standard now those are alternating. Honestly, I am not fully aware of the advantages, it came from the definition in treams, do not know if there are advantages for you in your way of storing.


> As I see, in your parameter sweep you plan to save in the same file only wavelength variation, and you do not have any method parameters that change for different wavelengths, so this is all in order.  But the sweep is somehow inner, not the outer dimension, this can cause problems. 


In computation, now the approach is to have computation/files/ group that includes the source files needed to repeat the computation, and computation/method_parameters group for all the additional parameters that are used by the software.  


The name of the material is easier to search in the database if specified as an attribute, so the group itself can be named however you like, but the attribute "name" which has the value "Ag" or "Silver"  or better both, would be the place to search for the material name. The reference to permittivity values(as a link to refractiveindex.info etc.) can be in the attribute "reference" of the group, could add also the interpolation method name in the attribute "interpolation" (those are some new aspects that came after discussion in the workshop).  Similarly in geometry, the group name should match the group name in materials/, and have the attribute  "shape" to set the name of the shape as you have it, datasets like radius etc are part of the group. You also specify shape as dataset with list of letters, this we did not have before, I do not know whether this brings extra advantages.  You are specifying in the name of the shape  that the spheroid is oblate or prolate, which we did not do before, but this is searchable under spheroid so is fine. "unit" can be added as an attribute to this group, if you do not provide a mesh. 


We have discussed the field for analytical_zeros with others, and for some time thought to maybe add it to geometry, but computation as it is now actually looks better. However, you store it as a flattened array of indices, and qp then refers to the modes? I suspect this was just for your personal use, sorry that I try to bring this to some standard that we document, but we probably would set some more explanatory naming convention. 

 

The issue with the mesh representation of the objects, I think, applies to you in the least sense, because you do not have any mesh discretization, however, we need to see for the database whether it is always required to have some geometry representation to avoid ambiguity in the position and orientation, for example if there is an arrangement of rotated spheroids to store. By default, no additional rotation/translation is assumed.


Final small remark, turned out the identifier for the database will be assigned after submission to the database, so uuid is not a requirement anymore, but could be good for local use, of course.


So far about the file, sorry for some new rearrangements, these will be solidified in the final document and repository.
Now regarding the wavelength sweep, I think there is no preference, and it's solely based on what you find reasonable. It is assumed that one would be able to search for a specific range of wavelengths in the database and obtain just that part for downloading, I do not know how well this will work in the end. 


Finally,  we were considering the possibility of adding a one-page subsection in the Generation section on how to create T-matrices using EBCM/SMARTIES and would like to inquire about your opinion on this matter.


