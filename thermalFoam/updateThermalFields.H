/*---------------------------------------------------------------------------*\
    updateThermalFields.H - Simplified Post-Erosion Update
    
    After mesh topology changes, this simply:
    1. Corrects boundary conditions on T
    2. Calls unified constitutive update
    3. Preserves boundary patch information
\*---------------------------------------------------------------------------*/

if (meshChanged)
{
    Info<< "Updating thermal fields on modified mesh..." << endl;
    
    // CRITICAL: Fix boundary conditions FIRST, before computing constitutive relations
    T.correctBoundaryConditions();
    theta.correctBoundaryConditions();
    
    // Update ALL constitutive relations in one unified call
    // This file is included, so all variables (Kth, CthEff, etc.) are available
    #include "updateConstitutive.H"
    
    // Boundary preservation logic for eroded patches
    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchI];
        const word& patchName = patch.name();
        
        if (patchName == erodedPatchName)
        {
            fvPatchScalarField& Tpatch = T.boundaryFieldRef()[patchI];
            if (Tpatch.type() == "codedMixed") continue;
            
            const labelList& faceCells = patch.faceCells();
            forAll(faceCells, faceI)
            {
                const label cellI = faceCells[faceI];
                if (originalPatchID[cellI] >= 0)
                {
                    label origPatchID = label(originalPatchID[cellI]);
                    if (origPatchID < mesh.boundaryMesh().size())
                    {
                        const word& origPatchName = mesh.boundaryMesh()[origPatchID].name();
                        if (erosionBoundaries.found(origPatchName))
                        {
                            Tpatch[faceI] = T[cellI];
                        }
                    }
                }
            }
        }
    }
    
    // Final boundary correction
    T.correctBoundaryConditions();
    
    Info<< "Thermal field update completed." << endl;
}
