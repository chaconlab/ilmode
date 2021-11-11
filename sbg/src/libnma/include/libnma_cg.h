#ifndef LIBNMA_CG_H_
#define LIBNMA_CG_H_

#include "nma.h"

// Creates a CA-model with first NH and last CO pseudo-atoms of each segment.
// Warning, atoms are not copied, they're just pointers to the original atoms.
// setmodel = true --> Sets masses in Occupancies and Number of electrons in Bfactors, otherwise left unchanged.
// equalmass = true --> Sets 1.0 masses to all CAs excepting NH and CO at segment endings
//                      (CA, NH and CO of the first and last segment residues will weight 0.5)
// equalmass = false --> Total residue mass is applied to all CAs excepting those at segment endings which will weight
//                       the Residue_mass-(NH_mass or CO_mass). Terminal NH and CO will weight 14 and 28, respectively.
// equalnelec = true --> Sets 1 electron to all CAs excepting NH and CO at segment endings
//                       (CA, NH and CO of the first and last segment residues will account for 0.5 electrons)
// equalnelec = false --> Total residue number of electrons is applied to all CAs excepting those at segment endings. The latter
//                        will have the total minus the N or CO number of electrons (7 or 14, respectively).
// setnelec = true --> Sets Number-of-electrons in Bfactors, otherwise they are left unchanged (for iModfit)
// setnelec = false --> Sets Number-of-electrons in Bfactors, otherwise they are left unchanged (for iMod's programs)
Macromolecule *cg_CA( Macromolecule *molNCAC, bool setmodel = true, bool equalmass = false, bool equalnelec = false, bool setnelec = false );


// Selects a CA-model with first NH and last CO pseudo-atoms of each segment.
// Warning, atoms are not copied, they're just pointers to the original atoms.
Macromolecule *select_cg_CA( Macromolecule *molNCAC);

// Sets masses to a NH,CA,CO-model (3-atoms/residue) like in CA-model (EXPERIMENTAL)
// equalmass = true --> sets 1.0 masses to all CAs excepting NH and CO at segment endings
//                   (unit mass is divided between NH or CO and their CAs)
// equalmass = false --> whole residue mass applied to all CAs excepting NH and CO at segment endings
//                   ( 15, 28 and Residue_mass-(NH_mass or CO_mass) for NH, CO and its CAs)
// allatoms = true --> ALL NH and CO will weight 15 and 28, respectively.
//                   ( CA's will weight the remaining residue mass )
// allatoms = false --> ALL allatoms-segment NH and CO will have 0.000001 mass, just to force spring linking.
//                   ( CA's will weight the residue mass )
// equalnelec = true --> Each atom will account for just one electron.
// equalnelec = false --> NH and CO pseudo-atoms will have the number of electrons of their heavy atoms.
//                        CA pseudo-atom will account for the corresponding remaining number of electrons.
// setnelec = true --> Sets Number-of-electrons in Bfactors, otherwise they are left unchanged (for iModfit)
// setnelec = false --> Sets Number-of-electrons in Bfactors, otherwise they are left unchanged (for iMod's programs)
void mass_NCAC( Macromolecule *mol, bool equalmass = false, bool allatoms = false, bool equalnelec = false, bool setnelec = false );

// CREATES a 3BB2R reduced model
//     Each residue is represented by 3 backbone atoms and 2 side-chain atoms.
//     Thus, 3 dihedral angles are needed for each residue (phi, psi, chi).
//     There are a few exceptions: for Ala, Gly and Pro,
//     and for the 1st and last residues.
// equalmass = true --> all masses set to 1.0 (default=false)
// equalnelec = true --> all charges set to 1.0 (default=false)
// setnelec = true --> Sets Number-of-electrons in Bfactors, otherwise they are left unchanged (for iModfit)
// setnelec = false --> Sets Number-of-electrons in Bfactors, otherwise they are left unchanged (for iMod's programs)
void cg_3BBR2( Macromolecule *mol, bool equalmass = false, bool equalnelec = false, bool setnelec = false);

//// CREATES a 3BB2R reduced model
////     Each residue is represented by 3 backbone atoms and 2 side-chain atoms.
////     Thus, 3 dihedral angles are needed for each residue (phi, psi, chi).
////     There are a few exceptions: for Ala, Gly and Pro,
////     and for the 1st and last residues.
//void cg_3BBR2( Macromolecule *mol );

// Sets masses to a Full-Atom-model
// equalmass = true --> all masses set to 1.0 (default=false)
// equalmass = false --> each atom will weight its own mass
// equalnelec = true --> all charges set to 1.0 (default=false)
// equalnelec = false --> each atom will account for its own number of electrons (default=false)
// setnelec = true --> Sets Number-of-electrons in Bfactors, otherwise they are left unchanged (for iModfit)
// setnelec = false --> Sets Number-of-electrons in Bfactors, otherwise they are left unchanged (for iMod's programs)
void mass_FA(Macromolecule *mol, bool equalmass = false, bool equalnelec = false, bool setnelec = false);

// Sets masses to a CA-only-model
// nomass = true --> masses = 1.0
// nomass = false --> each atom will weight its mass
void mass_CAonly(Macromolecule *mol, bool nomass = false);

// Places the "O" atom in the proper sp2-hybridation place
// (from a 3BB2R-model) (for visualization purposes only!)
void backboner(Macromolecule *mol);

#endif /*LIBNMA_CG_H_*/
