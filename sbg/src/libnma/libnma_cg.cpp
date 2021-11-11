/************************************************************************
*                     LIBRARY: libnma_cg                                *
*************************************************************************
* Program is part of the ADP package URL: http://sbg.cib.csic.es        *
* (c) Jose Ramon Lopez-Blanco (Mon), Jose Ignacio Garzon and            *
* Pablo Chacon. CSIC-CIB's Structural Bioinformatics Group (2004-2010)  *
*************************************************************************
*                                                                       *
*   Coarse-Graining models library.                                     *
*   (It takes into account Muliple-Chains and different CG-models)      *
*                                                                       *
*************************************************************************
* This library is free software; you can redistribute it and/or modify  *
* it under the terms of the GNU General Public License as published by  *
* the Free Software Foundation; either version 2 of the License, or     *
* (at your option) any later version.                                   *
************************************************************************/

// ######################################################################################
//
// Coarse-Graining NMA related library
//
// ######################################################################################
//#include <nma.h>
#include <libnma_cg.h>

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
Macromolecule *cg_CA( Macromolecule *molNCAC, bool setmodel, bool equalmass, bool equalnelec, bool setnelec )
{
	bool debug = false;
	pdbIter *iter,*iter2,*iterA;
	float imass; // current mass
	float inelec; // current number of electrons
	int num_res;
	Residue *res;
	Atom *at;

	if( setmodel ) // Setting Masses and Number of electrons according to the CA-only model
	{
		iter = new pdbIter( molNCAC ); // iter to screen fragments (residues)
		if(debug)
			if(equalmass)
				printf( "Msg(cg_CA): CA's will weight 1.0\n");
			else
				printf( "Msg(cg_CA): CA's will weight their residue mass.\n");

		for ( iter->pos_segment = 0; !iter->gend_segment(); iter->next_segment() )
		{
			iter2 = new pdbIter( iter->get_segment() );
			num_res = iter2->num_fragment();
			for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() ) // iter frags
			{
				res = ( Residue * ) iter2->get_fragment();
				iterA = new pdbIter( res ); // screens residue atoms

				// Applying whole residue mass to CA-atom
				iterA->pos_atom = 1; // CA's position
				at = iterA->get_atom();
				// Obtaining current masses and number of electrons
				if(iter2->pos_fragment == 0)
				{
					if(equalmass)
						imass = 0.5; // the 1.0 mass is divided between CA and NH atoms
					else
						imass = AA[ resnum_from_resname( res->getName() ) ].mass - 14; // whole residue mass - N

					if(setnelec)
					{
						if(equalnelec)
							inelec = 0.5; // one electron is divided between CA and NH atoms
						else
							inelec = AA[ resnum_from_resname( res->getName() ) ].nelec - 7; // whole residue number of electrons minus N's ones
//						fprintf(stderr,"setnelec enabled, inelec= %f\n",inelec);
					}
				}
				else if(iter2->pos_fragment == num_res-1)
				{
					if(equalmass)
						imass = 0.5; // the 1.0 mass is divided between CA and CO
					else
						imass = AA[ resnum_from_resname( res->getName() ) ].mass - 28; // whole residue mass - CO

					if(setnelec)
					{
						if(equalnelec)
							inelec = 0.5; // one electron is divided between CA and CO atoms
						else
							inelec = AA[ resnum_from_resname( res->getName() ) ].nelec - 14; // whole residue number of electrons minus CO's ones
					}
				}
				else
				{
					if(equalmass)
						imass = 1.0; // forcing unit masses
					else
						imass = AA[ resnum_from_resname( res->getName() ) ].mass; // whole residue mass - CO

					if(setnelec)
					{
						if(equalnelec)
							inelec = 1.0; // one electron
						else
							inelec = AA[ resnum_from_resname( res->getName() ) ].nelec; // whole residue number of electrons minus CO's ones
//						fprintf(stderr,"setnelec enabled, inelec= %f\n",inelec);
					}
				}
				at->setPdbocc( imass ); // mass
				if(setnelec)
					at->setPdbfact( inelec ); // number of electrons (for pdb to map translation)

				// Applying zero mass to N-atom
				iterA->pos_atom = 0; // N's position (virtual)
				at = iterA->get_atom();
				if(iter2->pos_fragment == 0)
				{
					if(equalmass)
						imass = 0.5; // the 1.0 mass is divided between CA and NH atoms
					else
						imass = 14; // N's mass

					if(setnelec)
					{
						if(equalnelec)
							inelec = 0.5; // one electron is divided between CA and NH atoms
						else
							inelec = 7; // N's number of electrons
					}
				}
				else
				{
					imass = 0.0;
					inelec = 0.0; // N's number of electrons
				}
				at->setPdbocc( imass ); // mass
				if(setnelec)
					at->setPdbfact( inelec ); // number of electrons (for pdb to map translation)

				// Applying zero mass to C-atom
				iterA->pos_atom = 2; // C's position (virtual)
				at = iterA->get_atom();
				if(iter2->pos_fragment == num_res-1)
				{
					if(equalmass)
						imass = 0.5; // the 1.0 mass is divided between CA and CO atoms
					else
						imass = 28; // CO's mass

					if(setnelec)
					{
						if(equalnelec)
							inelec = 0.5; // one electron is divided between CA and CO atoms
						else
							inelec = 14; // CO's number of electrons
					}
				}
				else
				{
					imass = 0.0;
					inelec = 0.0; // N's number of electrons
				}
				at->setPdbocc( imass ); // mass
				if(setnelec)
					at->setPdbfact( inelec ); // number of electrons (for pdb to map translation)

				delete iterA;
			}
			delete iter2;
		}
		delete iter;
	}

	// This creates a "by-reference" CA-model copy from the initial Macromolecule with (First-NH + CA's + Last-CO) model
	Macromolecule *mol;
	mol = new Macromolecule();
	Protein *pro,*pro2;
	Chain *ch,*ch2;
	Segment *seg,*seg2;
	Residue *frag,*frag2;
	pdbIter *iter_ch;
	pdbIter *iter_seg;
	pdbIter *iter_frag;
	pdbIter *iter_at;
	pdbIter *iter_prot = new pdbIter( molNCAC );
	for( iter_prot->pos_molecule = 0; !iter_prot->gend_molecule(); iter_prot->next_molecule() )
	{
		pro = (Protein *) iter_prot->get_molecule();
		pro2 = new Protein( pro->getName(), pro->getIdNumber() );
		iter_ch = new pdbIter( pro ); // iters current protein
		for( iter_ch->pos_chain = 0; !iter_ch->gend_chain(); iter_ch->next_chain() )
		{
			ch = iter_ch->get_chain();
			ch2 = new Chain( ch->getName(), ch->getIdNumber() );
			iter_seg = new pdbIter( ch ); // iters current chain
			for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
			{
				seg = iter_seg->get_segment();
				seg2 = new Segment( seg->getName(), seg->getIdNumber() );
				iter_frag = new pdbIter( seg ); // iters current segment
				num_res = iter_frag->num_fragment();
				for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
				{
					frag = (Residue *) iter_frag->get_fragment();
					frag2 = new Residue( frag->getName(), frag->getIdNumber(), frag->get_pos(), frag->get_letter() );
					iter_at = new pdbIter( frag );
					if(iter_frag->pos_fragment == 0) // if first segment fragment
					{
						if(debug)
							printf("%d %d %d %d  Adding NH + CA\n",iter_prot->pos_molecule
									,iter_ch->pos_chain,iter_seg->pos_segment	,iter_frag->pos_fragment);
						// get NH
						iter_at->pos_atom = 0; // NH
						frag2->add(iter_at->get_atom());
						// get CA
						iter_at->pos_atom = 1; // CA
						frag2->add(iter_at->get_atom());
					}
					else if(iter_frag->pos_fragment == num_res-1) // if last segment fragment
					{
						if(debug)
							printf("%d %d %d %d  Adding CA + CO\n",iter_prot->pos_molecule
									,iter_ch->pos_chain,iter_seg->pos_segment	,iter_frag->pos_fragment);
						// get CA
						iter_at->pos_atom = 1; // CA
						frag2->add(iter_at->get_atom());
						// get CO
						iter_at->pos_atom = 2; // CO
						frag2->add(iter_at->get_atom());
					}
					else
					{
						if(debug)
							printf("%d %d %d %d  Adding CA\n",iter_prot->pos_molecule
									,iter_ch->pos_chain,iter_seg->pos_segment,iter_frag->pos_fragment);
						// get CA
						iter_at->pos_atom = 1; // CA
						frag2->add(iter_at->get_atom());
					}
					if(debug)
						printf("TMOL = %d\n",frag->getMolType());
					delete iter_at;
					seg2->add(frag2);
				}
				delete iter_frag;
				ch2->add(seg2);
			}
			delete iter_seg;
			pro2->add(ch2);
		}
		mol->add(pro2);
		delete iter_ch;
	}
	delete iter_prot;
	return mol;
}

// Selects a CA-model with first NH and last CO pseudo-atoms of each segment.
// Warning, atoms are not copied, they're just pointers to the original atoms.
Macromolecule *select_cg_CA( Macromolecule *molNCAC)
{
	bool debug = false;
	int num_res,i;

	// SELECTING ALL (by reference, "atoms" will not be copied, only referenced)
	// Conditions
	Condition * all = new Condition( -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 );
	Conditions * all2 = new Conditions();
	all2->add( all );

	// First-NH + CA's + Last-CO model
	Macromolecule *mol;
//	mol = new Macromolecule(molNCAC); // copy initial molecule
	mol = molNCAC->select(all2,true); // select by reference
//	mol->writePDB("all.pdb");

	pdbIter *iter_seg = new pdbIter(mol);
	pdbIter *iter_frag;
	Residue *frag;

	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		iter_frag = new pdbIter( iter_seg->get_segment() ); // iters current segment
		num_res = iter_frag->num_fragment();
		for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			frag = (Residue *)iter_frag->get_fragment();
			if(iter_frag->pos_fragment == 0) // if first segment fragment
				for(i = frag->getLimit()-1; i > 1; i--) // delete all not (NH or CA)
					frag->remove(i);
			else if(iter_frag->pos_fragment == num_res-1) // if last segment fragment
			{
				for(i = frag->getLimit()-1; i>2 ; i--) // delete all above (NH,CA,C)
					frag->remove(i);
				frag->remove(0);
			}
			else // middle segment
			{
				for(i = frag->getLimit()-1; i>1 ; i--) // delete all above (NH,CA)
					frag->remove(i);
				frag->remove(0); // delete NH
			}
		}
		delete iter_frag;
	}
	delete iter_seg;
	return mol;
}

// Selects a CA-model with first NH and last CO pseudo-atoms of each segment.
// Warning, atoms are not copied, they're just pointers to the original atoms.
Macromolecule *select_cg_CA_old( Macromolecule *molNCAC)
{
	bool debug = false;
	int num_res;

	// First-NH + CA's + Last-CO model
	Macromolecule *mol;
	mol = new Macromolecule();
	Protein *pro,*pro2;
	Chain *ch,*ch2;
	Segment *seg,*seg2;
	Residue *frag,*frag2;
	pdbIter *iter_ch;
	pdbIter *iter_seg;
	pdbIter *iter_frag;
	pdbIter *iter_at;
	pdbIter *iter_prot = new pdbIter( molNCAC );
	for( iter_prot->pos_molecule = 0; !iter_prot->gend_molecule(); iter_prot->next_molecule() )
	{
		pro = (Protein *) iter_prot->get_molecule();
		pro2 = new Protein( pro->getName(), pro->getIdNumber() );
		iter_ch = new pdbIter( pro ); // iters current protein
		for( iter_ch->pos_chain = 0; !iter_ch->gend_chain(); iter_ch->next_chain() )
		{
			ch = iter_ch->get_chain();
			ch2 = new Chain( ch->getName(), ch->getIdNumber() );
			iter_seg = new pdbIter( ch ); // iters current chain
			for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
			{
				seg = iter_seg->get_segment();
				seg2 = new Segment( seg->getName(), seg->getIdNumber() );
				iter_frag = new pdbIter( seg ); // iters current segment
				num_res = iter_frag->num_fragment();
				for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
				{
					frag = (Residue *) iter_frag->get_fragment();
					frag2 = new Residue( frag->getName(), frag->getIdNumber(), frag->get_pos(), frag->get_letter() );
					iter_at = new pdbIter( frag );
					if(iter_frag->pos_fragment == 0) // if first segment fragment
					{
						if(debug)
							printf("%d %d %d %d  Adding NH + CA\n",iter_prot->pos_molecule
									,iter_ch->pos_chain,iter_seg->pos_segment	,iter_frag->pos_fragment);
						// get NH
						iter_at->pos_atom = 0; // NH
						frag2->add(iter_at->get_atom());
						// get CA
						iter_at->pos_atom = 1; // CA
						frag2->add(iter_at->get_atom());
					}
					else if(iter_frag->pos_fragment == num_res-1) // if last segment fragment
					{
						if(debug)
							printf("%d %d %d %d  Adding CA + CO\n",iter_prot->pos_molecule
									,iter_ch->pos_chain,iter_seg->pos_segment	,iter_frag->pos_fragment);
						// get CA
						iter_at->pos_atom = 1; // CA
						frag2->add(iter_at->get_atom());
						// get CO
						iter_at->pos_atom = 2; // CO
						frag2->add(iter_at->get_atom());
					}
					else
					{
						if(debug)
							printf("%d %d %d %d  Adding CA\n",iter_prot->pos_molecule
									,iter_ch->pos_chain,iter_seg->pos_segment,iter_frag->pos_fragment);
						// get CA
						iter_at->pos_atom = 1; // CA
						frag2->add(iter_at->get_atom());
					}
					if(debug)
						printf("TMOL = %d\n",frag->getMolType());
					delete iter_at;
					seg2->add(frag2);
				}
				delete iter_frag;
				ch2->add(seg2);
			}
			delete iter_seg;
			pro2->add(ch2);
		}
		mol->add(pro2);
		delete iter_ch;
	}
	delete iter_prot;
	return mol;
}

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
void mass_NCAC( Macromolecule *mol, bool equalmass, bool allatoms, bool equalnelec, bool setnelec )
{
	bool debug = false;
	// Setting Masses according to the CA-only model
	pdbIter *iterA,*iter,*iter2;
	int num_res;
	iter = new pdbIter( mol ); // iter to screen fragments (residues)
	Residue *res;
	Atom *at;

	float imass; // current mass
	float inelec; // current number of electrons

	if(debug)
		if(equalmass) // mass = 1.0
			if(allatoms)
				printf( "Msg(mass_NCAC): CA, NH and CO, will weight 1.0\n");
			else
				printf( "Msg(mass_NCAC): CA mass = 1.0, NH = 0.0, and CO = 0.0\n");
		else // account for residue mass
			if(allatoms)
				printf( "Msg(mass_NCAC): NH and CO pseudo-atoms will weight 14 and 28, and CA the remaining residue mass.\n");
			else
				printf( "Msg(mass_NCAC): CA's will weight whole residue mass, NH = 0.0, and CO = 0.0\n");

	for ( iter->pos_segment = 0; !iter->gend_segment(); iter->next_segment() )
	{
		iter2 = new pdbIter( iter->get_segment() );
		num_res = iter2->num_fragment();
		for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() ) // iter frags
		{
			res = ( Residue * ) iter2->get_fragment();
			iterA = new pdbIter( res ); // screens residue atoms

			// Applying masses to CA
			iterA->pos_atom = 1; // CA's position
			at = iterA->get_atom();

			if(equalmass)
				imass = 1.0; // atomic mass
			else
			{
				if(allatoms)
					imass = AA[ resnum_from_resname( res->getName() ) ].mass - 43; // CA-pending mass
				else if(iter2->pos_fragment == 0) // if first residue
					imass = AA[ resnum_from_resname( res->getName() ) ].mass - 14; // CA-pending mass
				else if(iter2->pos_fragment == num_res-1) // if first residue
					imass = AA[ resnum_from_resname( res->getName() ) ].mass - 28; // CA-pending mass
				else
					imass = AA[ resnum_from_resname( res->getName() ) ].mass; // CA= full residue mass
			}
			at->setPdbocc( imass ); // mass

			if(setnelec)
			{
				if(equalnelec)
					inelec = 1.0; // number of electrons
				else
				{
					if(allatoms)
						inelec = AA[ resnum_from_resname( res->getName() ) ].nelec - 21; // CA-pending number of electrons
					else if(iter2->pos_fragment == 0) // if first residue
						inelec = AA[ resnum_from_resname( res->getName() ) ].nelec - 7; // CA-pending number of electrons
					else if(iter2->pos_fragment == num_res-1) // if first residue
						inelec = AA[ resnum_from_resname( res->getName() ) ].nelec - 14; // CA-pending number of electrons
					else
						inelec = AA[ resnum_from_resname( res->getName() ) ].nelec; // total number of electrons
				}
				at->setPdbfact( inelec ); // number of electrons
			}

			// Applying mass to N-atom
			iterA->pos_atom = 0; // N's position (virtual)
			at = iterA->get_atom();
			if(equalmass)
				imass = 1.0; // atomic mass
			else
				if(iter2->pos_fragment == 0 || allatoms) // if first residue or "allatoms" masses enabled
					imass = 14.0; // NH mass
				else
					imass = 0.000001; // != 0 forces spring connection
			at->setPdbocc( imass );

			if(setnelec)
			{
				if(equalnelec)
					inelec = 1.0; // atomic mass
				else
					if(iter2->pos_fragment == 0 || allatoms) // if first residue or "allatoms" masses enabled
						inelec = 7; // N number of electrons
					else
						inelec = 0.0; // no electrons
				at->setPdbfact( inelec );
			}

			// Applying mass to C-atom
			iterA->pos_atom = 2; // C's position (virtual)
			at = iterA->get_atom();
			if(equalmass)
				imass = 1.0; // atomic mass
			else
				if(iter2->pos_fragment == num_res-1 || allatoms) // if last residue
					imass = 28; // CO mass
				else
					imass = 0.000001; // != 0 forces spring connection
			at->setPdbocc( imass );

			if(setnelec)
			{
				if(equalnelec)
					inelec = 1.0; // atomic mass
				else
					if(iter2->pos_fragment == num_res-1 || allatoms) // if last residue
						inelec = 14; // CO number of electrons
					else
						inelec = 0.0; // no electrons
				at->setPdbfact( inelec );
			}

			delete iterA;
		}
		delete iter2;
	}
	delete iter;
}

// CREATES a 3BB2R reduced model
//     Each residue is represented by 3 backbone atoms and 2 side-chain atoms.
//     Thus, 3 dihedral angles are needed for each residue (phi, psi, chi).
//     There are a few exceptions: for Ala, Gly and Pro,
//     and for the 1st and last residues.
// equalmass = true --> all masses set to 1.0 (default=false)
// equalnelec = true --> all charges set to 1.0 (default=false)
// setnelec = true --> Sets Number-of-electrons in Bfactors, otherwise they are left unchanged (for iModfit)
// setnelec = false --> Sets Number-of-electrons in Bfactors, otherwise they are left unchanged (for iMod's programs)
void cg_3BBR2( Macromolecule *mol, bool equalmass, bool equalnelec, bool setnelec)
{
	Residue * res;
	int j, resn, atoms, atoms_R, atom_res;
	char * at_name;
	Segment * seg;
	Atom * atom;
	Tcoor pos;
	pdbIter * iter1, * iter2;
	float cx, cy, cz; // center of masses
	int num_res=0;
	float masstotal, mass, massCH, massCH2,massCH3, massNH, massCO;
	float nelectotal, charge, nelecCH, nelecCH2, nelecCH3, nelecNH, nelecCO;

	// Electronic Charges (more EM-like!)
	if(setnelec)
		if(equalnelec)
		{
			nelecCH = 1.0;
			nelecCH2 = 1.0;
			nelecCH3 = 1.0;
			nelecNH = 1.0;
			nelecCO = 1.0;
		}
		else
		{
			nelecCH = 6 + 1;
			nelecCH2 = nelecCH + 1;
			nelecCH3 = nelecCH2 + 1;
			nelecNH = 7 + 1;
			nelecCO = 6 + 8;
		}

	// Masses
	if(equalmass)
	{
		massCH = 1.0;
		massCH2 = 1.0;
		massCH3 = 1.0;
		massNH = 1.0;
		massCO = 1.0;
	}
	else
	{
		massCH = 12.011 + 1.00797;
		massCH2 = massCH + 1.00797;
		massCH3 = massCH2 + 1.00797;
		massNH = 14.00674 + 1.00797;
		massCO = 12.011 + 15.9994;
	}

	//Iterador para recorrer segmentos
	iter1 = new pdbIter( mol );

	//Bucle para recorrer segmentos
	for ( iter1->pos_segment = 0; !iter1->gend_segment(); iter1->next_segment() )
	{
		seg = ( Segment * ) iter1->get_segment();
		iter2 = new pdbIter( seg );

		//Bucle para recorrer residuos
		for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() )
		{
			res = ( Residue * ) iter2->get_fragment();
			resn = resnum_from_resname( res->getName() );
			num_res++;
			//Iterador sobre los atomos del residuo
			pdbIter * iter3 = new pdbIter( res );

			switch (resn)
			{
			case GLY:
				if (iter3->num_atom() < 3)
				{
					printf( " Residue  %s%d only %d atoms\n", res->getName(),num_res, iter3->num_atom());
					exit(1);
				}
				iter3->pos_atom = 0; //NH
				(( Atom * ) iter3->get_atom())->setPdbocc(massNH);
				if(setnelec)
					(( Atom * ) iter3->get_atom())->setPdbfact(nelecNH);

				iter3->next_atom(); // CA+H
				(( Atom * ) iter3->get_atom())->setPdbocc(massCH);
				if(setnelec)
					(( Atom * ) iter3->get_atom())->setPdbfact(nelecCH);

				iter3->next_atom(); // CO
				(( Atom * ) iter3->get_atom())->setPdbocc(massCO);
				if(setnelec)
					(( Atom * ) iter3->get_atom())->setPdbfact(nelecCO);

				// Deleting remaining un-necessary atoms
				res->remove(3);

				break;
			case ALA:
				if (iter3->num_atom() < 4)
				{
					printf( " Residue  %s%d only %d atoms\n", res->getName(),num_res, iter3->num_atom());
					exit(1);
				}
				iter3->pos_atom = 0; //NH
				(( Atom * ) iter3->get_atom())->setPdbocc(massNH);
				if(setnelec)
					(( Atom * ) iter3->get_atom())->setPdbfact(nelecNH);

				iter3->next_atom(); // CA+HA
				(( Atom * ) iter3->get_atom())->setPdbocc(massCH);
				if(setnelec)
					(( Atom * ) iter3->get_atom())->setPdbfact(nelecCH);

				iter3->next_atom(); // CO
				(( Atom * ) iter3->get_atom())->setPdbocc(massCO);
				if(setnelec)
					(( Atom * ) iter3->get_atom())->setPdbfact(nelecCO);

				//				iter3->next_atom(); // CB
				res->exchange( 3, 4 );   // exchange O for CB

				iter3->pos_atom = 4; // O in "res" structure, but CB in "iter3"
				iter3->get_atom()->setPdbocc(massCH3);
				if(setnelec)
					iter3->get_atom()->setPdbfact(nelecCH3);

				// Deleting remaining un-necessary atoms
				res->remove(4); // in "res" this is O

				break;

			default:
				if (iter3->num_atom() < 5)
				{
					printf( " Residue  %s%d only %d atoms\n", res->getName(),num_res, iter3->num_atom());
					exit(1);
				}

				iter3->pos_atom = 0; //NH
				(( Atom * ) iter3->get_atom())->setPdbocc(massNH);
				if(setnelec)
					(( Atom * ) iter3->get_atom())->setPdbfact(nelecNH);

				iter3->next_atom(); // CA+HA
				(( Atom * ) iter3->get_atom())->setPdbocc(massCH);
				if(setnelec)
					(( Atom * ) iter3->get_atom())->setPdbfact(nelecCH);

				iter3->next_atom(); // CO
				(( Atom * ) iter3->get_atom())->setPdbocc(massCO);
				if(setnelec)
					(( Atom * ) iter3->get_atom())->setPdbfact(nelecCO);

				//				iter3->next_atom(); // CB
				res->exchange( 3, 4 );   // exchange O for CB

				iter3->pos_atom = 4; // O in "res" structure, but CB in "iter3"
				if (( resn == ILE ) ||( resn == THR ) ||( resn == VAL ))
				{
					(( Atom * ) iter3->get_atom())->setPdbocc(massCH); // Setting CH mass for CB pseudo-atom (pos 4)
					if(setnelec)
						(( Atom * ) iter3->get_atom())->setPdbfact(nelecCH); // Setting CH number of electrons for CB pseudo-atom (pos 4)
				}
				else
				{
					(( Atom * ) iter3->get_atom())->setPdbocc(massCH2); // Setting CH2 mass for CB pseudo-atom (pos 4)
					if(setnelec)
						(( Atom * ) iter3->get_atom())->setPdbfact(nelecCH2); // Setting CH2 number of electrons for CB pseudo-atom (pos 4)
				}

				// 5th R2 --> center of mass of the rest of the chain..
				cx = 0;
				cy = 0;
				cz = 0;
				atoms = 0;

				masstotal=0;
				nelectotal=0;
				atom_res=5;
				atoms_R=0;
				for ( iter3->pos_atom = 5; !iter3->gend_atom(); iter3->next_atom() )
				{
					atom = ( Atom * ) iter3->get_atom();
					at_name = atom->getName();
					atom->getPosition(pos);
					atoms_R++;
					//printf( " %d ++%s++ %f %f %f\n", iter3->pos_atom, at_name, pos[0],pos[1],pos[2]);
					if ( strncmp( at_name, " H  ",4 ) )
						if ( strncmp( at_name, " HN ",4 ) )
							if ( strncmp( at_name, " HA ",4 ) )
								if ( strncmp( at_name, " HB ",4 ) )
									if ( strncmp( at_name, "1HB ",4 ) )
										if ( strncmp( at_name, "2HB ",4 ) )
										{ // whether it is not those Hydrogens...
											if(equalmass)
												mass= 1.0;
											else
												mass= ( atom->getElement() )->weight; // atomic weight
											if(equalnelec)
												charge = 1.0;
											else
												charge = ( atom->getElement() )->number; // atomic number (charge)
											masstotal+=mass;
											nelectotal+=charge;
											cx += pos[0]*mass;
											cy += pos[1]*mass;
											cz += pos[2]*mass;
											atoms++;
											//printf (" %s Mass %f pos  %f %f %f\n",at_name,  mass, pos[0],pos[1],pos[2]);
											//getchar();
										}
				}
				pos[0] = cx / masstotal;
				pos[1] = cy / masstotal;
				pos[2] = cz / masstotal;
				//printf (" res %s masstotal %f atomos R %d\n",res->getName(), masstotal, atoms_R);

				// 5th rest of the chain..
				iter3->pos_atom = 3; // "O" ???
				atom = ( Atom * ) iter3->get_atom();
				atom->setPosition(pos);
				atom->setPdbocc(masstotal);
				if(setnelec)
					atom->setPdbfact(nelectotal);
				iter3->~pdbIter();
				//printf( " mass %f cxyz %f %f %f atoms %d %d\n",
				//masstotal, pos[0],pos[1],pos[2],atoms, atoms_R);
				if (AA[resn].nheavyatoms!=atom_res+atoms_R)
					printf(" Warning res %s have %d atoms instead of %d\n",res->getName(), atom_res+atoms_R, AA[resn].nheavyatoms);
				for ( j = 0 ; j < atoms_R; j++)
					res->remove(5);
				break;
			}
		}
		delete iter2;
	}
	delete iter1;
}

// Sets masses to a Full-Atom-model
// equalmass = true --> all masses set to 1.0 (default=false)
// equalmass = false --> each atom will weight its own mass
// equalnelec = true --> all charges set to 1.0 (default=false)
// equalnelec = false --> each atom will account for its own number of electrons (default=false)
// setnelec = true --> Sets Number-of-electrons in Bfactors, otherwise they are left unchanged (for iModfit)
// setnelec = false --> Sets Number-of-electrons in Bfactors, otherwise they are left unchanged (for iMod's programs)
void mass_FA(Macromolecule *mol, bool equalmass, bool equalnelec, bool setnelec)
{
	bool debug = false;
	pdbIter *iter = new pdbIter( mol ); // iter to screen atoms
	Atom *at;
	float imass,inelec;

	if(debug)
		if(equalmass)
			printf( "Msg(mass_FA): Full-Atom model --> Atoms will weight 1.0\n");
		else
			printf( "Msg(mass_FA): Full-Atom model --> Atoms will weight their atomic mass.\n");

	for ( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() )
	{
		at = ( Atom * ) iter->get_atom();
		if(equalmass)
			imass = 1.0; // forcing unit masses
		else
			imass = at->getElement()->weight; // atomic mass
		at->setPdbocc( imass ); // mass

		if(setnelec)
		{
			if(equalnelec)
				inelec = 1.0; // one electron
			else
				inelec = ( at->getElement() )->number; // number of electrons
			at->setPdbfact( inelec ); // number of electrons
		}
		if(debug)
			printf("atom %d  id= %s  mass= %f\n",iter->pos_atom,at->getName(),imass);
	}

	delete iter;
}

// Sets masses to a CA-only-model
// equalmass = true --> masses = 1.0
// equalmass = false --> each CA will weight its residue mass
void mass_CAonly(Macromolecule *mol, bool equalmass)
{
	bool debug = false;
	pdbIter *iter = new pdbIter( mol ); // iter to screen atoms
	pdbIter *iter2 = new pdbIter( mol ); // iter to screen fragments
	Residue *res;
	Atom *at;
	float imass;

	if(debug)
		if(equalmass)
			printf( "Msg(mass_CAonly): CA-only model --> Atoms will weight 1.0\n");
		else
			printf( "Msg(mass_CAonly): CA-only model --> Atoms will weight their whole redidue mass.\n");

	for ( iter2->pos_fragment = 0, iter->pos_atom = 0; !iter2->gend_fragment();
			iter2->next_fragment(), iter->next_atom() )
	{
		res = ( Residue * ) iter2->get_fragment();
		if(equalmass)
			imass = 1.0; // forcing unit masses
		else
			imass = AA[ resnum_from_resname( res->getName() ) ].mass; // Residue mass
		at = ( Atom * ) iter->get_atom();
		at->setPdbocc( imass ); // mass
		at->setPdbfact( imass ); // number of electrons (DO IT BETTER LATER!!!)
	}
	delete iter;
	delete iter2;
}

// Places the "O" atom in the proper sp2-hybridation place
// (from a 3BB2R-model) (for visualization purposes only!)
void backboner(Macromolecule *mol)
{
	pdbIter *it_seg,*it_res,*it_at;
	Segment *seg;
	Residue *res;
	Atom *at,*atO;
	int resn,num_res;
	Tcoor posN,posCA,posC,posO,e1,e2;
	float module;

	it_seg = new pdbIter( mol );
	for( it_seg->pos_segment = 0; !it_seg->gend_segment(); it_seg->next_segment() )
	{
		seg = ( Segment * ) it_seg->get_segment();
		it_res = new pdbIter( seg );
		num_res = it_res->num_fragment(); // gets segment's number of residues
		for ( it_res->pos_fragment = 0; !it_res->gend_fragment(); it_res->next_fragment() )
		{
			res = ( Residue * ) it_res->get_fragment();
			resn = resnum_from_resname( res->getName() );
			it_at = new pdbIter( res );

			it_at->pos_atom = 0; //NH
			at = it_at->get_atom();
			at->getPosition(posN);

			// Re-placing "O" atom
			if(it_res->pos_fragment > 0 && it_res->pos_fragment < num_res-1 ) // not first or last residue
			{
				// e1 = vector C-CA
				e1[0] = posCA[0] - posC[0];
				e1[1] = posCA[1] - posC[1];
				e1[2] = posCA[2] - posC[2];
				// module=1   vector C-CA
				module = sqrt( e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2] );
				e1[0] /= module;
				e1[1] /= module;
				e1[2] /= module;

				// e2 = vector C-N
				e2[0] = posN[0] - posC[0];
				e2[1] = posN[1] - posC[1];
				e2[2] = posN[2] - posC[2];
				// module=1   vector C-N
				module = sqrt( e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2] );
				e2[0] /= module;
				e2[1] /= module;
				e2[2] /= module;

				// Now e1' = e1 + e2
				e1[0] += e2[0];
				e1[1] += e2[1];
				e1[2] += e2[2];
				// module=1   vector C-CA ^ C-N (bisectriz-vector)
				module = sqrt( e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2] );
				e1[0] /= module;
				e1[1] /= module;
				e1[2] /= module;

				// C=O canonical bond-length = 1.23 A
				posO[0] = posC[0] -1.23 * e1[0];
				posO[1] = posC[1] -1.23 * e1[1];
				posO[2] = posC[2] -1.23 * e1[2];

				atO->setPosition( posO ); // setting last O-atom position
				atO->setPdbocc(16.0);
				atO->setPdbfact(8.0);
			}

			if( it_res->pos_fragment == num_res-1 ) // if last residue within segment
			{
				// O-atoms in the last segment residue... are different!
				// Under construction...........
			}

			it_at->next_atom(); // CA+HA
			at = it_at->get_atom();
			at->getPosition(posCA);

			it_at->next_atom(); // CO
			at = it_at->get_atom();
			at->getPosition(posC);
			at->setPdbocc(12.011); // sets only C-atom mass

			// if GLY or ALA, an extra O-atom must be added!
			if( resn == GLY || resn == ALA )
			{
			      atO = new Atom( Table_Elements::getElement( "O " ), posN, 0.0, " O   ", 1 , 15.9994, 8.0 );
			      res->add( atO ); // adding new atom
			}
			else // standard residue...
			{
				it_at->pos_atom = 4; // O
				atO = it_at->get_atom();
			}
			// O-atom position will be updated later...
			// ...when the next residue N-atom position is known!

			delete it_at;
		}
		delete it_res;
	}
	delete it_seg;
}
