# Notes

- codedMixed fails setting the value at the boundary of the patch. Could be
an issue of U = fvc::grad(h) or a problem with how codedMixed update the 
values on the patch. It is easier to daisychain cases that flip between a fixedValue and
hotstart modifying h to a fixedGradient. 