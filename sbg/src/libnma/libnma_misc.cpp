/*************************************************************************
 *                     LIBRARY: libnma_hessian                           *
 *************************************************************************
 * Program is part of the ADP package URL: http://sbg.cib.csic.es        *
 * (c) Jose Ramon Lopez-Blanco (Mon), Jose Ignacio Garzon and            *
 * Pablo Chacon. CSIC-CIB's Structural Bioinformatics Group (2004-2010)  *
 *************************************************************************
 *                                                                       *
 *   Miscelanea NMA related rutines library, to:                         *
 *   	-Make interacting pair of atoms (ipas)                           *
 *      -Create residue properties (proper)                              *
 *      -Translate modes and fixation arrays into different CG-models    *
 *      -Create segment properties (addrot)                              *
 *      -Etc...                                                          *
 *   (It takes into account Muliple-Chains and different CG-models)      *
 *                                                                       *
 *************************************************************************
 * This library is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 ************************************************************************/

#include <libtools/include/Memstuff.h>

#include "libnma_io.h"
#include "libnma_misc.h"
CRandomMersenne *rg; // Mersenne Twister global object

#define DEVEL

// Creates contacts list (IPA) directly from a Macromolecule ("mol")
void make_ipas_new(Macromolecule *mol, twid **p_decint, int *p_nipa, float cutoff)
{
  bool debug = false;
  int nipa, index, k, l;
  double d;
  twid *decint;
  decint = *p_decint;
  nipa = 0; // counts the number of interacting pseudo-atom pairs
pdbIter *iterA,*iterB;
iterA = new pdbIter(mol);
iterB = new pdbIter(mol);
Atom *atA,*atB;
Tcoor rA,rB;

  if(debug)
	  printf("Msg(make_ipas_new): Creating Interacting Pairs of pseudo-Atoms (IPAs) list.\n");

  for(iterA->pos_atom=0; !iterA->gend_atom(); iterA->next_atom() )
  {
	  atA = iterA->get_atom();
		  atA->getPosition(rA);
		  for(iterB->pos_atom=iterA->pos_atom+1; !iterB->gend_atom(); iterB->next_atom() )
		  {
			  atB = iterB->get_atom();
				  atB->getPosition(rB);
				  d = sqrt( pow(rA[0]-rB[0],2) + pow(rA[1]-rB[1],2) + pow(rA[2]-rB[2],2) );

				  if(d <= cutoff)
				  {
					  nipa++; // Counts number of Interacting Pairs of Atoms
					  if( !(decint = ( twid * ) realloc( decint, nipa * sizeof( twid ) ) ) ) // resizes contact list-structure
					  {
						  printf("Sorry, \"decint\" memmory re-allocation failed!\n");
						  exit(1);
					  }
					  decint[nipa - 1].k = iterA->pos_atom; // k-pseudo-atom index (i-atom index)
					  decint[nipa - 1].l = iterB->pos_atom; // l-pseudo-atom index (j-atom index)
					  decint[nipa - 1].d = d; // set distance
					  decint[nipa - 1].C = 0.0; // force constant will be set in the future
				  }
		  }
  }
  delete iterA;
  delete iterB;

  if(debug)
	  printf( "Msg(make_ipas_new): Number of Interacting Pairs of pseudo-Atoms (NIPAs): %d\n", nipa );

  *p_nipa=nipa; // outputs "nipas"
  *p_decint=decint;
}

// Creates contacts list (IPA) directly from a Macromolecule ("mol")
// (The "masses", i.e. occupancies, with zero mass will not be accounted for.)
void make_ipas_new0(Macromolecule *mol, twid **p_decint, int *p_nipa, float cutoff)
{
  bool debug = false;
  int nipa, index, k, l;
  double d;
  twid *decint;
  decint = *p_decint;
  nipa = 0; // counts the number of interacting pseudo-atom pairs
pdbIter *iterA,*iterB;
iterA = new pdbIter(mol);
iterB = new pdbIter(mol);
Atom *atA,*atB;
Tcoor rA,rB;

  if(debug)
	  printf("Msg(make_ipas_new): Creating Interacting Pairs of pseudo-Atoms (IPAs) list.\n");

  for(iterA->pos_atom=0; !iterA->gend_atom(); iterA->next_atom() )
  {
	  atA = iterA->get_atom();
	  if(atA->getPdbocc() != 0.0) // if it's not a virtual atom (modified 24/11/2009)
	  {
		  atA->getPosition(rA);
		  for(iterB->pos_atom=iterA->pos_atom+1; !iterB->gend_atom(); iterB->next_atom() )
		  {
			  atB = iterB->get_atom();
			  if(atB->getPdbocc() != 0.0)  // if it's not a virtual atom
			  {
				  atB->getPosition(rB);
				  d = sqrt( pow(rA[0]-rB[0],2) + pow(rA[1]-rB[1],2) + pow(rA[2]-rB[2],2) );

				  if(d <= cutoff)
				  {
					  nipa++; // Counts number of Interacting Pairs of Atoms
					  decint = ( twid * ) realloc( decint, nipa * sizeof( twid ) ); // resizes contact list-structure
					  decint[nipa - 1].k = iterA->pos_atom; // k-pseudo-atom index (i-atom index)
					  decint[nipa - 1].l = iterB->pos_atom; // l-pseudo-atom index (j-atom index)
					  decint[nipa - 1].d = d; // set distance
					  decint[nipa - 1].C = 0.0; // force constant will be set in the future
				  }
			  }
		  }
	  }
  }
  delete iterA;
  delete iterB;

  if(debug)
	  printf( "Msg(make_ipas_new): Number of Interacting Pairs of pseudo-Atoms (NIPAs): %d\n", nipa );

  *p_nipa=nipa; // outputs "nipas"
  *p_decint=decint;
}


// Multiplies by "factor" every ipa force constant ("C") if both atoms belong to different molecules.
void modify_intermolec_ipas(Macromolecule *mol, twid *decint, int nipa, float factor)
{
  bool debug = true;
  pdbIter *iterA,*iterB;
  iterA = new pdbIter(mol);
  iterB = new pdbIter(mol);

  Molecule *molecA,*molecB;

  if(debug)
	  printf("Msg(modify_ipas): Modifying Interacting Pairs of pseudo-Atoms (IPAs) list.\n");

  for(int i=0; i<nipa;i++)
  {
	  iterA->pos_atom = decint->k;
	  iterB->pos_atom = decint->l;
	  molecA = ((Molecule *)((Chain *)((Segment *)((Fragment *)(iterA->get_atom())->getFather())->getFather())->getFather())->getFather());
	  molecB = ((Molecule *)((Chain *)((Segment *)((Fragment *)(iterB->get_atom())->getFather())->getFather())->getFather())->getFather());

	  printf("Msg(modify_ipas): Interaction atoms: %5d %5d  molec: %2d %2d multiplied by %f\n",decint->k,decint->l,molecA->getIdNumber(),molecB->getIdNumber(),factor);
	  if(molecA->getIdNumber() != molecB->getIdNumber())
	  {
		  if(debug)
			  printf("Msg(modify_ipas): Interaction atoms: %5d %5d  molec: %2d %2d multiplied by %f\n",decint->k,decint->l,molecA->getIdNumber(),molecB->getIdNumber(),factor);
		  decint->C *= factor;
	  }
  }
  delete iterA;
  delete iterB;

  if(debug)
	  printf( "Msg(make_ipas_new): Number of Interacting Pairs of pseudo-Atoms (NIPAs): %d\n", nipa );

}


// Creates contacts list (IPA) directly from a Macromolecule ("mol")
// Contacts according to a predefined set of functions (depending on Topology and Secondary Structure)
void make_ipasTS(Macromolecule *mol, twid **p_decint, int *p_nipa, float cutoff, TSfunc *funcs, int nfuncs, char *ss)
{
	bool debug = false;
	int nipa, index, k, l, top, a_res=0, b_res=0, a_atom=0, b_atom=0;
	double d,fc;
	twid *decint;
	decint = *p_decint;

	if(debug)
		printf("Msg(make_ipasTS): Creating Interacting Pairs of pseudo-Atoms (IPAs) list.\n");

	nipa = 0; // counts the number of interacting pseudo-atom pairs
	pdbIter *iterA,*iterB,*iterAres,*iterBres,*iterAseg,*iterBseg;
	iterAseg = new pdbIter(mol); // to iterate at segment level
	iterBseg = new pdbIter(mol);
	Tcoor rA,rB;
	bool contact = false;
	bool null_ss = false;

	if(ss==NULL) // Just topology... ( 'X' vs. 'X' )
		null_ss = true;

	for( iterAseg->pos_segment=0; !iterAseg->gend_segment(); iterAseg->next_segment() )
	{
		iterAres = new pdbIter( iterAseg->get_segment() );
		for( iterAres->pos_fragment=0; !iterAres->gend_fragment(); iterAres->next_fragment() )
		{
			iterA = new pdbIter( iterAres->get_fragment() ); // iter to screen atoms inside a residue
			for(iterA->pos_atom = 0; !iterA->gend_atom(); iterA->next_atom() )
			{
				(iterA->get_atom())->getPosition(rA);

				b_res = a_res;
				b_atom = a_atom+1;
				for( iterBseg->pos_segment=iterAseg->pos_segment; !iterBseg->gend_segment(); iterBseg->next_segment() )
				{
					iterBres = new pdbIter( iterBseg->get_segment() );

					if( iterBseg->pos_segment == iterAseg->pos_segment )
						iterBres->pos_fragment = iterAres->pos_fragment;
					else
						iterBres->pos_fragment = 0;

					for( ; !iterBres->gend_fragment(); iterBres->next_fragment() )
					{
						iterB = new pdbIter( iterBres->get_fragment() ); // iter to screen atoms inside a residue

						if( iterBres->pos_fragment == iterAres->pos_fragment && iterBseg->pos_segment == iterAseg->pos_segment )
							iterB->pos_atom = iterA->pos_atom+1; // current A-atom with the B-atoms from the same B-residue as the A-residue which A-atom belongs to.
						else
							iterB->pos_atom = 0; // current A-atom with the remaining B-atoms in B-residues different from current A-residue.

						for( ; !iterB->gend_atom(); iterB->next_atom() )
						{
							(iterB->get_atom())->getPosition(rB);
							d = sqrt( pow(rA[0]-rB[0],2) + pow(rA[1]-rB[1],2) + pow(rA[2]-rB[2],2) );

							if(d <= cutoff)
							{
								top = iterBres->pos_fragment - iterAres->pos_fragment; // topology level

								for(int n=0; n<nfuncs && !contact && !null_ss; n++)
								{
									// SS matching
									if(ss[a_res] == funcs[n].i && ss[b_res] == funcs[n].j)
									{
										if(funcs[n].t < 0) // non-bonding matching
										{
											fc = inv_exp(funcs[n].a,d,funcs[n].b,funcs[n].c); // setting Force Constants
											contact = true; // k-l pair contacted
										}
										else if(top == funcs[n].t) // Topology matching
										{
											fc = inv_exp(funcs[n].a,d,funcs[n].b,funcs[n].c); // setting Force Constants
											contact = true; // k-l pair contacted
										}
									}
								}

								// If k-l pair doesn't belong to the SS-matching functions (i.e. the XX functions)
								for(int n=0; n<nfuncs && !contact; n++)
								{
									// XX matching
									if( funcs[n].i == 'X' && funcs[n].j == 'X' )
									{
										if(funcs[n].t < 0) // non-bonding matching
										{
											fc = inv_exp(funcs[n].a,d,funcs[n].b,funcs[n].c); // setting Force Constants
											contact = true; // k-l pair contacted
										}
										// Topology matching
										else if(top == funcs[n].t)
										{
											fc = inv_exp(funcs[n].a,d,funcs[n].b,funcs[n].c); // setting Force Constants
											contact = true; // k-l pair contacted
										}
									}
								}

								if(contact)
								{
									decint = ( twid * ) realloc( decint, (nipa+1) * sizeof( twid ) ); // resizes contact list-structure
									decint[nipa].k = a_atom;; // k-pseudo-atom index
									decint[nipa].l = b_atom; // l-pseudo-atom index
									decint[nipa].d = d; // set distance
									decint[nipa].C = fc; // set force constant
									// debug info
									if(debug)
										printf("resA= %4d (%c) atomA= %4d  resB= %4d (%c) atomB= %4d --> C= %8.5f\n",a_res,ss[a_res],a_atom,b_res,ss[b_res],b_atom,decint[nipa].C);
									nipa++; // Counts number of Interacting Pairs of Atoms
								}
								contact = false;
							}
							b_atom++; // secuential atom index
						}
						b_res++; // B-residue index
						delete iterB;
					}
					delete iterBres;
				}
				a_atom++; // secuential atom index
			}
			a_res++; // A-residue index
			delete iterA;
		}
		delete iterAres;
	}
	delete iterAseg;
	delete iterBseg;

	if(debug)
		printf( "Msg(make_ipasTS): Number of Interacting Pairs of pseudo-Atoms (NIPAs): %d\n", nipa );

	*p_nipa=nipa; // outputs "nipas"
	*p_decint=decint;
}

//	LAS REGLAS DEL MIX:
//
//	1) Sea un par de residuos i , j:
//	S=|i-j| (distancia secuencial)
//	Si S<=3 , entonces Kij=Cseq / S� donde Cseq=60 kcal/mol.A�
//	Si S>3, entonces Kij=(Ccart/dij)^6, donde Ccart=6 "
//
//	2) Ademas hay un cutoff igual a:
//	Si N=length > 50, cutoff=8 Angtroms
//	Si N>50, entonces
//	cutoff=integer(Sl*log(N)-In)
//	donde Sl=3 , In=2.8.
//
// "MIXED" contacting method. (WARNING: only fully valid for single chain and CA models)
// Creates contacts list (IPA) and assigns force constants given a Macromolecule ("mol")
// (The "masses", i.e. occupancies, with zero mass will not be accounted for.)
void make_ipas_mix0(Macromolecule *mol, twid **p_decint, int *p_nipa)
{
	int nipa, index, dist_s, num_res,i,j,offset;
	double d,cutoff;
	twid *decint;
	decint = *p_decint;
	nipa = 0; // counts the number of interacting pseudo-atom pairs
	pdbIter *iterA,*iterB,*iterA_frag,*iterB_frag;
	iterA_frag = new pdbIter(mol);
	iterB_frag = new pdbIter(mol);
	Atom *atA,*atB;
	Tcoor rA,rB;
	Fragment *fragA,*fragB;
	double Cseq = 60.0; // kcal/mol.A^2
	double Ccart = 6.0;
	double power = 6.0; // exponential term
	double Sl = 3.0;
	double In = 2.8;

	// Setting cutoff: cutoff=integer(Sl*log(N)-In) (log--> Natural log)
	num_res = mol->get_num_fragments();
	if(num_res<50)
		cutoff = 8;
	else
		cutoff = (int) ( Sl * log((float)num_res) - In );
	//fprintf(stderr,"Msg(make_ipas_mix0): Contacting by \"mixed\" method: cutoff= %f\n",cutoff);

	i=0;
	for(iterA_frag->pos_fragment = 0; !iterA_frag->gend_fragment(); iterA_frag->next_fragment() )
	{
		fragA = iterA_frag->get_fragment();
		iterA = new pdbIter(fragA); // iterate "A" fragment atoms

		for(iterA->pos_atom=0; !iterA->gend_atom(); iterA->next_atom() )
		{
			atA = iterA->get_atom();
			if(atA->getPdbocc() != 0.0) // if it's not a virtual atom (modified 24/11/2009)
			{
				atA->getPosition(rA);

				j=i;
				for(iterB_frag->pos_fragment = iterA_frag->pos_fragment;
						!iterB_frag->gend_fragment(); iterB_frag->next_fragment() )
				{
					fragB = iterB_frag->get_fragment();
					iterB = new pdbIter(fragB); // iterate "B" fragment atoms
					dist_s = iterB_frag->pos_fragment - iterA_frag->pos_fragment; // Secuential distance (fragment level)
//				    fprintf(stderr,"fragA= %d  fragB= %d  S= %d\n",iterA_frag->pos_fragment,iterB_frag->pos_fragment,dist_s);
					if(iterB_frag->pos_fragment == iterA_frag->pos_fragment)
						offset=iterA->pos_atom;
					else
						offset=0;
					for(iterB->pos_atom=offset; !iterB->gend_atom(); iterB->next_atom() )
					{
						atB = iterB->get_atom();
						if(atA != atB) // if atoms are different!
						{
							if(atB->getPdbocc() != 0.0)  // if it's not a virtual atom
							{
								atB->getPosition(rB);
								d = sqrt( pow(rA[0]-rB[0],2) + pow(rA[1]-rB[1],2) + pow(rA[2]-rB[2],2) );
// 							    fprintf(stderr,"\tatomA= %d  atomB= %d  i= %d  j= %d  d= %f\n",iterA->pos_atom,iterB->pos_atom,i,j,d);

								if(dist_s==0) // Needed by non CA-only models
								{	//	If S=0 , then Kij=Cseq  , where Cseq=60 kcal/mol.A^2 (MON added)
									nipa++;
									decint = ( twid * ) realloc( decint, nipa * sizeof( twid ) ); // resizes contact list-structure
									decint[nipa - 1].k = i; // k-pseudo-atom index (i-atom index)
									decint[nipa - 1].l = j; // l-pseudo-atom index (j-atom index)
									decint[nipa - 1].d = d; // set distance
									decint[nipa - 1].C = (double) Cseq;
								}
								else if(dist_s<=3)
								{	//	If S<=3 , then Kij=Cseq / S�  , where Cseq=60 kcal/mol.A^2
									nipa++;
									decint = ( twid * ) realloc( decint, nipa * sizeof( twid ) ); // resizes contact list-structure
									decint[nipa - 1].k = i; // k-pseudo-atom index (i-atom index)
									decint[nipa - 1].l = j; // l-pseudo-atom index (j-atom index)
									decint[nipa - 1].d = d; // set distance
									decint[nipa - 1].C = (double) Cseq/(dist_s*dist_s); // force constant will be set in the future
								}
								else if(d <= cutoff)
								{	//	If S>3, then Kij=(Ccart/dij)^6  , where Ccart=6
									nipa++;
									decint = ( twid * ) realloc( decint, nipa * sizeof( twid ) ); // resizes contact list-structure
									decint[nipa - 1].k = i; // k-pseudo-atom index (i-atom index)
									decint[nipa - 1].l = j; // l-pseudo-atom index (j-atom index)
									decint[nipa - 1].d = d; // set distance
									decint[nipa - 1].C = pow(Ccart/d,power); // force constant will be set in the future
								}
							}
						}
						j++;
					}
					delete iterB;
				}
			}
			i++;
		}
		delete iterA;
	}
	delete iterA_frag;
	delete iterB_frag;

	*p_nipa=nipa; // outputs "nipas"
	*p_decint=decint; // and "ipas"
}

// "MIXED" contacting method. (WARNING: only fully valid for single chain and CA models)
// Fills a contact_matrix with force constants given a Macromolecule ("mol")
void cont_matrix_mix(Macromolecule *mol, double *cont_matrix, int num_atoms, int *p_nipa)
{
	int nipa, index, dist_s, num_res,i,j,offset;
	double d,cutoff;
	nipa = 0; // counts the number of interacting pseudo-atom pairs
	pdbIter *iterA,*iterB,*iterA_frag,*iterB_frag;
	iterA_frag = new pdbIter(mol);
	iterB_frag = new pdbIter(mol);
	Atom *atA,*atB;
	Tcoor rA,rB;
	Fragment *fragA,*fragB;
	double Cseq = 60.0; // kcal/mol.A^2
	double Ccart = 6.0;
	double power = 6.0; // exponential term
	double Sl = 3.0;
	double In = 2.8;

	// Setting cutoff: cutoff=integer(Sl*log(N)-In) (log--> Natural log)
	num_res = mol->get_num_fragments();
	if(num_res<50)
		cutoff = 8;
	else
		cutoff = (int) ( Sl * log((float)num_res) - In );
	//fprintf(stderr,"Msg(make_ipas_mix0): Contacting by \"mixed\" method: cutoff= %f\n",cutoff);

	i=0;
	for(iterA_frag->pos_fragment = 0; !iterA_frag->gend_fragment(); iterA_frag->next_fragment() )
	{
		fragA = iterA_frag->get_fragment();
		iterA = new pdbIter(fragA); // iterate "A" fragment atoms

		for(iterA->pos_atom=0; !iterA->gend_atom(); iterA->next_atom() )
		{
			atA = iterA->get_atom();
			if(atA->getPdbocc() != 0.0) // if it's not a virtual atom (modified 24/11/2009)
			{
				atA->getPosition(rA);

				j=i;
				for(iterB_frag->pos_fragment = iterA_frag->pos_fragment;
						!iterB_frag->gend_fragment(); iterB_frag->next_fragment() )
				{
					fragB = iterB_frag->get_fragment();
					iterB = new pdbIter(fragB); // iterate "B" fragment atoms
					dist_s = iterB_frag->pos_fragment - iterA_frag->pos_fragment; // Secuential distance (fragment level)
//				    fprintf(stderr,"fragA= %d  fragB= %d  S= %d\n",iterA_frag->pos_fragment,iterB_frag->pos_fragment,dist_s);
					if(iterB_frag->pos_fragment == iterA_frag->pos_fragment)
						offset=iterA->pos_atom;
					else
						offset=0;
					for(iterB->pos_atom=offset; !iterB->gend_atom(); iterB->next_atom() )
					{
						atB = iterB->get_atom();
						if(atA != atB) // if atoms are different!
						{
							if(atB->getPdbocc() != 0.0)  // if it's not a virtual atom
							{
								atB->getPosition(rB);
								d = sqrt( pow(rA[0]-rB[0],2) + pow(rA[1]-rB[1],2) + pow(rA[2]-rB[2],2) );

								if(dist_s==0) // Needed by non CA-only models
								{	//	If S=0 , then Kij=Cseq  , where Cseq=60 kcal/mol.A^2 (MON added)
									cont_matrix[ num_atoms * i + j - 1 - i * ( i + 3 ) / 2 ] = (double) Cseq;
									nipa++;
								}
								else if(dist_s<=3)
								{	//	If S<=3 , then Kij=Cseq / S�  , where Cseq=60 kcal/mol.A^2
									cont_matrix[ num_atoms * i + j - 1 - i * ( i + 3 ) / 2 ] = (double) Cseq/(dist_s*dist_s);
									nipa++;
								}
								else if(d <= cutoff)
								{	//	If S>3, then Kij=(Ccart/dij)^6  , where Ccart=6
									cont_matrix[ num_atoms * i + j - 1 - i * ( i + 3 ) / 2 ] = pow(Ccart/d,power);
									nipa++;
								}
							}
						}
						j++;
					}
					delete iterB;
				}
			}
			i++;
		}
		delete iterA;
	}
	delete iterA_frag;
	delete iterB_frag;

	*p_nipa=nipa; // outputs "nipas"
}

// Yang,...,Jernigan's "pfENM". PNAS (2009)
// Fills a contact_matrix with force constants given a Macromolecule ("mol")
void cont_matrix_pfENM(Macromolecule *mol, double *cont_matrix, int num_atoms, int *p_nipa)
{
	int nipa, index, num_res,i,j;
	double d,cutoff;
	nipa = 0; // counts the number of interacting pseudo-atom pairs
	pdbIter *iterA,*iterB;
	iterA = new pdbIter(mol);
	iterB = new pdbIter(mol);
	Atom *atA,*atB;
	Tcoor rA,rB;
	double k= 1.0; // kcal/mol.A^2
	double power = 2.0; // exponential term

	// Setting cutoff: cutoff=integer(Sl*log(N)-In) (log--> Natural log)
	num_res = mol->get_num_fragments();
	//fprintf(stderr,"Msg(cont_matrix_pfENM): Contacting by \"pfENM\" method.\n");

	for(i=iterA->pos_atom=0; !iterA->gend_atom(); iterA->next_atom(), i++ )
	{
		atA = iterA->get_atom();
		if(atA->getPdbocc() != 0.0) // if it's not a virtual atom (modified 24/11/2009)
		{
			atA->getPosition(rA);
			for(j=iterB->pos_atom=iterA->pos_atom+1; !iterB->gend_atom(); iterB->next_atom(), j++ )
			{
				atB = iterB->get_atom();
				if(atB->getPdbocc() != 0.0)  // if it's not a virtual atom
				{
					atB->getPosition(rB);
					d = sqrt( pow(rA[0]-rB[0],2) + pow(rA[1]-rB[1],2) + pow(rA[2]-rB[2],2) );

					cont_matrix[ num_atoms * i + j - 1 - i * ( i + 3 ) / 2 ] = k/pow(d,power);
					nipa++;
				}
			}
		}
	}
	delete iterA;
	delete iterB;

	*p_nipa=nipa; // outputs "nipas"
}

// Kovaks's "inverse exponential" (Rueda et al.)
// Fills a contact_matrix with force constants given a Macromolecule ("mol")
void cont_matrix_Kovacs(Macromolecule *mol, double *cont_matrix, int num_atoms, int *p_nipa)
{
	int nipa, index, num_res,i,j;
	double d,cutoff;
	nipa = 0; // counts the number of interacting pseudo-atom pairs
	pdbIter *iterA,*iterB;
	iterA = new pdbIter(mol);
	iterB = new pdbIter(mol);
	Atom *atA,*atB;
	Tcoor rA,rB;
	//  Kij= k((3.8/d_ij)^6)
	double k= 40.0; // Rueda et. al
	double d0=3.8;
	double power = 6.0; // exponential term

	// Setting cutoff: cutoff=integer(Sl*log(N)-In) (log--> Natural log)
	num_res = mol->get_num_fragments();
	//fprintf(stderr,"Msg(cont_matrix_pfENM): Contacting by \"pfENM\" method.\n");

	for(i=iterA->pos_atom=0; !iterA->gend_atom(); iterA->next_atom(), i++ )
	{
		atA = iterA->get_atom();
		if(atA->getPdbocc() != 0.0) // if it's not a virtual atom (modified 24/11/2009)
		{
			atA->getPosition(rA);
			for(j=iterB->pos_atom=iterA->pos_atom+1; !iterB->gend_atom(); iterB->next_atom(), j++ )
			{
				atB = iterB->get_atom();
				if(atB->getPdbocc() != 0.0)  // if it's not a virtual atom
				{
					atB->getPosition(rB);
					d = sqrt( pow(rA[0]-rB[0],2) + pow(rA[1]-rB[1],2) + pow(rA[2]-rB[2],2) );

					cont_matrix[ num_atoms * i + j - 1 - i * ( i + 3 ) / 2 ] = k*pow(d0/d,power);
					nipa++;
				}
			}
		}
	}
	delete iterA;
	delete iterB;

	*p_nipa=nipa; // outputs "nipas"
}

// Computes some residue properties (Full-Atom & Multi-Chain & RNA/DNA & SMOL)
// Creates an array to relate atoms and the rigid units which atoms belongs to (unat).
void properMFA(Macromolecule *molr,tri **p_props,int **p_unat,int type,int model)
{
	bool debug=false;
	int i,resn=0,num_res,num_atom;
	tri *props;
// MON: "unat" not used... REMOVE! (MEMORY MAY BE NEEDED ANYWAYS...)
	int *unat;
	pdbIter *iter_frag; // Iterator to screen fragments
	pdbIter *iter_seg; // Iterator to screen segments
	pdbIter *iter_at; // Iterator to screen atoms
	int index_res = 0;

	if(debug)
		printf("Msg(properMFA): Debugging properMFA!\n");

	iter_seg = new pdbIter( molr, true, true, true, true ); // Iterator to screen segments

	num_res = iter_seg->num_fragment();
	if( !(props = ( tri * ) malloc( num_res * sizeof( tri ) ) ) )
	{
		printf("Msg(proper): Unable to allocate memory!\n");
		exit(1);
	}

	num_atom = molr->get_num_atoms();
// MON: "unat" not used... REMOVE! (MEMORY MAY BE NEEDED ANYWAYS...)
	if( !(unat = ( int * ) malloc( num_atom * sizeof( int ) ) ) )
	{
		printf("Msg(proper): Unable to allocate memory!\n");
		exit(1);
	}
	for(i=0;i<num_atom;i++)
		unat[i] = 0; // initialization

	props[0].k1 = 0; // index of first atom of first residue


	Segment *seg;
	Residue *res;
	TMOL fragtype;

	//	Atom * atom;
	// Screening segments
	for ( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
//fprintf(stderr,"Breakpoint, at segment %d\n",iter_seg->pos_segment);
//exit(0);
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)
		if(debug)
			printf ("Processing segment %d (%d res):\n",iter_seg->pos_segment,num_res);

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

		// Initializing some fragment properties (# atoms, # dihedrals, etc...)
		for ( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() ) // screen residues
		{
			i=iter_frag->pos_fragment; // shorter
			res = ( Residue * ) iter_frag->get_fragment();
//			fragtype = res->getMolType(); // This can be outside,...
			if(fragtype != tmol_smol)
				resn = resnum_from_resname( res->getName() );
			else
				resn = 666;

			iter_at = new pdbIter(res); // iter residue atoms
			props[index_res].nat = iter_at->num_atom(); // fragment number of atoms

//			if( iter_frag->pos_fragment > 0 )  // MON: Bug that worked because "k1" was not generally used...
			if(index_res > 0)
				props[index_res].k1 = props[index_res-1].k1 + props[index_res-1].nat; // index of first atom in current residue

//			fprintf(stderr,"index_res= %4d  nat= %1d  nan= %1d  k1= %4d\n", index_res, props[index_res].nat, props[index_res].nan, props[index_res].k1);

//		    fprintf(stderr,"%6d resn: %d  fragtype= %d  ADE=%d\n", iter_frag->pos_fragment, resn,fragtype,ADE);

			if( fragtype == tmol_protein )
			{
				if(num_res == 1 && model == 1) // If single residue segment
					props[index_res].nan = 0; // Nothing "flexible" should be done here.
				else
				{
					if ( resn == ALA ) // No-Chi
						props[index_res].nan=2;   // # dihedral angles
					else if ( resn == GLY ) // No-Chi
						props[index_res].nan = 2;
					else if ( resn == PRO )
						props[index_res].nan = 1;
					else
					{
						if(type == 0) // WITHOUT CHIs (like if it were an ALA with 5 atoms...)
							props[index_res].nan=2;   // # dihedral angles
						else if (type == 2 || type == 1) // (Phi,Chi,Psi)
							props[index_res].nan = 3;
					}
					// Check if: first (non proline) or last...
					if ( ( iter_frag->pos_fragment == 0 && resn != PRO ) || (iter_frag->pos_fragment == num_res - 1 && model == 1) )
						props[index_res].nan -= 1;
				}
			} // END-protein
			// RNA/DNA
			else if( resn == GUA || resn == ADE || resn == CYT || resn == URA
					|| resn == DGUA || resn == DADE || resn == DCYT || resn == DTHY )
			{
//			    fprintf(stderr,"%6d resn: %d  fragtype= %d  Is a nucleotide\n", iter_frag->pos_fragment, resn,fragtype);
				if(i==num_res-1) // last segment fragment
					props[index_res].nan = 3;   // no-epsilon and no-eta
				else // not last segment fragment
					props[index_res].nan = 5;   // full dihedrals
				if( type == 2 ) // if CHI is present
					props[index_res].nan++;
			} // END-RNA/DNA
			// SMOL - Small MOLecules (ligands)
			else if( fragtype == tmol_smol )
			{
				props[index_res].nan = 0; // Nothing "flexible" should be done here.
			}
			else // NOT FOUND MOL-TYPE
			{
				printf("Msg(properMFA): Sorry, Molecular-type (TMOL) not identified!\nForcing exit!\n");
				exit(2);
			}

			if(debug)
				printf("frag= %d (%3s) nat= %d nangles= %d\n",index_res,(iter_frag->get_fragment())->getName(),props[index_res].nat,props[index_res].nan); // residue number of atoms

//			fprintf(stderr,"index_res= %4d  nat= %1d  nan= %1d  k1= %4d\n", index_res, props[index_res].nat, props[index_res].nan, props[index_res].k1);

			index_res++; // counts residues

			delete iter_at;
		}
		delete(iter_frag);
	}
	iter_seg->clean_virtual();
	delete(iter_seg);
	*p_props = props; // returns allocated struct
	*p_unat = unat;
}

// Computes some residue properties (CA-model & Multi-Chain)
// "p_props" stores the properities.
// "p_unat" stores the unit which an atom belongs to.
// *WARNING "unat" is not properly filled here! (you need the "fix" array!)
void properCA(Macromolecule *molr,tri **p_props, int **p_unat)
{
	int i,resn=0,num_res,num_atom;
	tri *props;
// MON: "u_at" and "n_un" not used... REMOVE!
	int n_at=0,n_un=0;
	int *unat;
	int index_res = 0;
	bool debug=false;

	num_res = molr->get_num_fragments();
	if( !(props = ( tri * ) malloc( num_res * sizeof( tri ) ) ) )
	{
		printf("Msg(proper): Unable to allocate memory!\n");
		exit(1);
	}

// MON: "unat" not used... REMOVE! (MEMORY MAY BE NEEDED ANYWAYS...)
	num_atom = molr->get_num_atoms();
	if( !(unat = ( int * ) malloc( num_atom * sizeof( int ) ) ) )
	{
		printf("Msg(proper): Unable to allocate memory!\n");
		exit(1);
	}
	for(i=0;i<num_atom;i++)
		unat[i] = 0; // initialization

//	// MON: newly added...
//	props[0].k1 = 0; // index of first atom of first residue

	// Screening segments
	Segment *seg;
	Residue *res;
	pdbIter *iter_frag; // Iterator to screen fragments
	pdbIter *iter_seg; // Iterator to screen segments
	iter_seg = new pdbIter( molr ); // Iterator to screen segments
	for ( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)
		if(debug)
			printf ("Msg(properCA): Processing segment %d (%d res):\n",iter_seg->pos_segment,num_res);

// MON: "n_un" not used... REMOVE!
		if( iter_seg->pos_segment != 0 )
			n_un ++; // the next segment breaks the last unit

		// Initializing some aa properties (# atoms, # dihedrals, defining units)
		for ( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() ) // screen residues
		{
			i=iter_frag->pos_fragment; // shorter
			res = ( Residue * ) iter_frag->get_fragment();
			resn = resnum_from_resname( res->getName() );

			props[index_res].nat = 1; //  atom always!
			// Note: NCAC-(3 atom)model does not need "nat"...

//			// MON: newly added...
//			if( iter_frag->pos_fragment > 0 )
//				props[index_res].k1 = props[index_res-1].k1 + props[index_res-1].nat; // index of first atom in current residue

			if ( resn != PRO ) // if not-PRO --> Phi & Psi
			{
				if(i==0 || num_res == 1)
				{
					props[index_res].nan = 0;
//					unat[ n_at ] = n_un;
				}
				else if(i==num_res-1)
				{
					props[index_res].nan = 0;
//					unat[ n_at ] = n_un;
				}
				else
				{
					props[index_res].nan = 2;
//					unat[ n_at ] = n_un;
				}
				if(debug)
					printf("Msg(properCA): res %4d Not-PRO: angles= %d\n",i,props[index_res].nan );
			}
			else // if PRO --> Psi only
			{
				if(i==0 || num_res == 1)
				{
					props[index_res].nan = 0;
//					unat[ n_at ] = n_un;
				}
				else if(i==num_res-1)
				{
					props[index_res].nan = 0;
//					unat[ n_at ] = n_un;
				}
				else
				{
					props[index_res].nan = 1;
//					unat[ n_at ] = n_un;
				}
				if(debug)
					printf("Msg(properCA): res %4d PRO: angles= %d\n",i,props[index_res].nan );
			}
// MON: "u_at" not used... REMOVE!
			n_at += props[index_res].nat; // accumulated number of atoms
//			n_un += props[index_res].nan; // accumulated number of unit
			index_res++; // counts residues
		}
		delete(iter_frag);
	}
	delete(iter_seg);
	*p_props = props; // returns allocated struct
	// MON: "unat" not used... REMOVE!
	*p_unat = unat;
}

// Creates two auxiliar arrays with segment properties (needed due to fixing):
//   addrot[#seg] --> true, if 3 additional rotations should be added due to fixing.
//   effseg[#seg] --> <int>, with the number of "effective segment" for #seg.
// (Allocates memory itself, if it's needed)
int seg_props(Macromolecule *mol, tri *props, bool *fix, int model, int type, bool **p_addrot, int **p_effseg)
{
	bool debug=false;
	pdbIter *iter_seg,*iter_frag;
	Segment *seg;
	Residue *res;
	int resn,num_seg,num_res,seg_atoms_old;
	int ieff=0,seg_atoms_current=0,res_index=0,mobile=0; // effective segment index
	int *effseg=NULL; // Effective segment indices
	bool *addrot=NULL; // Should be added ROTATIONs? (with fixing)
	bool may_have_rot = false;
	bool split=false;
	int j=0;
	TMOL fragtype;

	iter_seg = new pdbIter( mol, true, true, true, true ); // Iterator to screen segments
	num_seg = iter_seg->num_segment();
	addrot = (bool *)malloc(sizeof(bool) * num_seg);
	effseg = (int *)malloc(sizeof(int) * num_seg);
	if( !addrot || !effseg)
	{
		printf("seg_props(): Sorry, memory allocation failure!\n");
		exit(1);
	}
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)

		// This fixes a bug...
		fragtype = seg->getMolType();

		// Checking single atom segments... SMOL's
		seg_atoms_old = seg_atoms_current;
		seg_atoms_current = seg->num_atoms();
		if(seg_atoms_old > 1)
			may_have_rot = true;

		if(debug)
			printf("seg= %d  num_res= %d  seg_atoms= %d (old=%d) may_have_rot= %d\n",iter_seg->pos_segment,num_res,seg_atoms_current,seg_atoms_old,may_have_rot);

		if(iter_seg->pos_segment != 0) // non-first segment
		{
			// TRANSLATION (always)
			for(int axis=0; axis<3; axis++)
			{
				if( fix==NULL || fix[j] ) // Check whether any TRANs variable it's mobile
				{
					split = true;
					mobile++; // count mobile DoFs
				}
				if(debug)
					if(fix==NULL)
						printf("seg_props(Msg): seg= %d  TRA(%d)= %d\n",iter_seg->pos_segment,axis,1);
					else
						printf("seg_props(Msg): seg= %d  TRA(%d)= %d\n",iter_seg->pos_segment,axis,fix[j]);
				j++;
			}

			// ROTATION
			if( seg_atoms_current != 1 && (seg_atoms_old != 1 || may_have_rot) )
				for(int axis=0; axis<3; axis++)
				{
					if( fix==NULL || fix[j] ) // Check whether any ROTs variable it's mobile
					{
						split = true;
						mobile++; // count mobile DoFs
					}
					if(debug)
						if(fix==NULL)
							printf("seg_props(Msg): seg= %d  ROT(%d)= %d\n",iter_seg->pos_segment,axis,1);
						else
							printf("seg_props(Msg): seg= %d  ROT(%d)= %d\n",iter_seg->pos_segment,axis,fix[j]);
					j++; // der[] index
				}

			// split segment
			if( split )
			{
				ieff++;
				split = false;
			}
		}

		// Screen ALL fragments
		for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			res = ( Residue * ) iter_frag->get_fragment();
//			fragtype = res->getMolType();
			if(fragtype != tmol_smol)
				resn = resnum_from_resname( res->getName() );
			else
				resn = 666;

			// if Protein fragment
			if( fragtype == tmol_protein )
			{
				// PHI --> NOT FIRST, NOT PRO, NOT CA-model's last residue
				if( !(iter_frag->pos_fragment == 0 || resn == PRO
						|| (iter_frag->pos_fragment == num_res-1 && (model == 0 || model == 3) ) ) )
				{
					if( fix==NULL || fix[j])
						mobile++; // count mobile DoFs
					if(debug)
						if(fix==NULL)
							printf("seg_props(Msg): seg= %d  frag= %d  PHI(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,1);
						else
							printf("seg_props(Msg): seg= %d  frag= %d  PHI(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,fix[j]);
					j++; // der[] index
				}

				//  LATERAL CHAIN-->CHI (NOT CA-model)
				if( type == 2 && model != 0 && model != 3)
					if( props[res_index].nan==3 ||
							(props[res_index].nan==2 &&
									(iter_frag->pos_fragment==0 ||
											(iter_frag->pos_fragment==num_res-1 && model==1) ) ) )
					{
						if( fix==NULL || fix[j])
							mobile++; // count mobile DoFs
						if(debug)
							if(fix==NULL)
								printf("seg_props(Msg): seg= %d  frag= %d  CHI(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,1);
							else
								printf("seg_props(Msg): seg= %d  frag= %d  CHI(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,fix[j]);
						j++; // der[] index
					}

				// "PSI-bond" (in Full-Atom, every residue has PSI) (NOT first residue in CA-model)
				if( !(iter_frag->pos_fragment == num_res-1
						|| (iter_frag->pos_fragment == 0 && (model==0 || model==3) ) )
						|| model==2 ) // Full-Atom allways has PSI
				{
					if( fix==NULL || fix[j])
						mobile++; // count mobile DoFs
					if(debug)
						if(fix==NULL)
							printf("seg_props(Msg): seg= %d  frag= %d  PSI(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,1);
						else
							printf("seg_props(Msg): seg= %d  frag= %d  PSI(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,fix[j]);
					j++; // der[] index
				}
			}
			// if RNA/DNA fragment
			else if( resn == GUA || resn == ADE || resn == CYT || resn == URA ||
					resn == DGUA || resn == DADE || resn == DCYT || resn == DTHY )
			{
				// ALPHA (bond between P and O5*)
				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					mobile++; // count mobile DoFs
				if(debug)
					if(fix==NULL)
						printf("seg_props(Msg): seg= %d  frag= %d  ALPHA(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,1);
					else
						printf("seg_props(Msg): seg= %d  frag= %d  ALPHA(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,fix[j]);
				j++; // der[] index

				// BETA (bond between O5* and C5*)
				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					mobile++; // count mobile DoFs
				if(debug)
					if(fix==NULL)
						printf("seg_props(Msg): seg= %d  frag= %d  BETA(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,1);
					else
						printf("seg_props(Msg): seg= %d  frag= %d  BETA(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,fix[j]);
				j++; // der[] index

				// GAMMA (bond between C5* and C4*)
				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					mobile++; // count mobile DoFs
				if(debug)
					if(fix==NULL)
						printf("seg_props(Msg): seg= %d  frag= %d  GAMMA(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,1);
					else
						printf("seg_props(Msg): seg= %d  frag= %d  GAMMA(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,fix[j]);
				j++; // der[] index

				//  LATERAL CHAIN --> CHI
				if(type == 2)
				{
					if( fix == NULL || fix[ j ] ) // Check whether the variable it's fixed
						mobile++; // count mobile DoFs
					if(debug)
						if(fix==NULL)
							printf("seg_props(Msg): seg= %d  frag= %d  CHI(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,1);
						else
							printf("seg_props(Msg): seg= %d  frag= %d  CHI(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,fix[j]);
					j++; // der[] index
				}

				// NOT LAST NUCLEOTIDE (4-IC's)
				if ( iter_frag->pos_fragment != num_res-1 )
				{
					// EPSILON (bond between C3* and O3*)
					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
						mobile++; // count mobile DoFs
					if(debug)
						if(fix==NULL)
							printf("seg_props(Msg): seg= %d  frag= %d  EPSILON(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,1);
						else
							printf("seg_props(Msg): seg= %d  frag= %d  EPSILON(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,fix[j]);
					j++; // der[] index

					// ZETA (bond between O3* and next-P)
					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
						mobile++; // count mobile DoFs
					if(debug)
						if(fix==NULL)
							printf("seg_props(Msg): seg= %d  frag= %d  ZETA(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,1);
						else
							printf("seg_props(Msg): seg= %d  frag= %d  ZETA(%d)= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,j,fix[j]);
					j++; // der[] index
				}
			}
			res_index++;
		}
		delete iter_frag;

		// Building effective segments array
		effseg[ iter_seg->pos_segment ] = ieff;
		if(debug)
			printf("segment= %2d  -->  effective= %2d\n",iter_seg->pos_segment,effseg[iter_seg->pos_segment]);
	}
	if(debug)
		printf("seg_props(Msg): Mobile(before)= %d\n",mobile);

	addrot[0]=false; // first segment never has ROTs
	seg_atoms_current = 0;
	may_have_rot = false;
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		// Checking single atom segments... SMOL's
		seg_atoms_old = seg_atoms_current;
		seg_atoms_current = seg->num_atoms();

		if(iter_seg->pos_segment != 0)
			if(    	(effseg[iter_seg->pos_segment] != effseg[iter_seg->pos_segment-1]) &&
					( (seg_atoms_current > 1) ||
							( iter_seg->pos_segment < num_seg-1 && effseg[iter_seg->pos_segment] ==
									effseg[iter_seg->pos_segment+1]) )	)
			{
				if(may_have_rot)
				{
					if(seg_atoms_current == 1)
					{
						mobile += 3; // count mobile DoFs
						addrot[iter_seg->pos_segment] = true; // add ROT when effective segment changes
					}
					else
						addrot[iter_seg->pos_segment] = false; // add ROT when effective segment changes
				}
				else
					addrot[iter_seg->pos_segment] = false; // nothing
				may_have_rot = true; // first time, none ROT should be added!
			}
			else
				addrot[iter_seg->pos_segment] = false; // nothing

		if(seg_atoms_old > 1 ||
				(iter_seg->pos_segment < num_seg-1 &&
						effseg[iter_seg->pos_segment] == effseg[iter_seg->pos_segment+1] ))
			may_have_rot = true;

		if(debug)
			printf("segment= %2d  -->  addrot= %d\n",iter_seg->pos_segment,(int)addrot[iter_seg->pos_segment]);
	}

	// Output
	*p_addrot = addrot;
	*p_effseg = effseg;

	if(debug)
		printf("Final DoFs= %d\n",mobile);

	return mobile;
}

// Creates a random mask of fixed internal coordinates according to a probability "prob"
// (Given a random number [0:1), an IC will be selected if "random>prob")
// "fix" array must be already allocated for full- "size" elements!
int fixRand(bool *fix,int size, float prob)
{

	bool debug = false;
	int mobile = 0;
	int i,r,temp;
	int *indices=NULL; // NULL is needed by realloc!

	if(debug)
		printf("Fixing Random (all coordinates)\n");

	indices = (int *)malloc(sizeof(int)*size);
	check_pointer(indices,"fixRand(), Indices array");

	// Creating array with coordinate indices
	for(i=0; i<size; i++)
		indices[i] = i;

	if(debug)
	{
		printf("fixRand> Coordinates found: %d\nfixRand> ",size);
		for(i=0; i<size; i++)
			printf("%d ",indices[i]);
		printf("\n");
	}

	// Now shuffle
	for(i=0; i<size; i++)
	{
//		r = size * (float) rg->Random(); // [0:dh_index-1) Playing dice with Mersenne!
		r = size * rg->Random(); // [0:dh_index-1) Playing dice with Mersenne!
		// swaping indices between "i" and "r"
		temp = indices[r];
		indices[r] = indices[i];
		indices[i] = temp;
	}

	// Showing shuffle
	if(debug)
	{
		printf("fixRand> Current Shuffle: ");
		for(i=0; i<size; i++)
			printf("%d ",indices[i]);
		printf("\n");
	}

	// Once shuffled... the first "prob-percent" dihedral coordinates will be fixed
	for(i=0; i<size; i++)
	{
		if(i < (int)size*(float)prob )
			fix[ indices[i] ] = 0; // fixed
		else
		{
			fix[ indices[i] ] = 1; // mobile
			mobile++; // count mobile coordinates
		}
	}
	free(indices);

	// Showing mask
	if(debug)
	{
		printf("fixRand> Current Mask: ");
		for(int i=0; i<size; i++)
			printf("%d ",fix[i]);
		printf("\nfixRand> mobile= %d\n",mobile);
	}

	return mobile;
}

//// Creates a random mask of fixed internal coordinates according to a probability "prob"
//// (Given a random number [0:1), an IC will be selected if "random>prob")
//// "fix" array must be already allocated for full- "size" elements!
//int fixRand(bool *fix,int size, float prob)
//{
////	CRandomMersenne * rg; // Mersenne Twister global object
////	extern CRandomMersenne *rg; // Mersenne Twister global object
//
//	bool debug=false;
//	float random;
//	int cont=0;
//
//	if(debug)
//		printf("Fixing Random (Rot/Tra/DH), size= %d\n",size);
//
//	for(int i=0; i<size; i++)
//	{
//		random = (float) rg->Random(); // [0:1) Playing dice with Mersenne!
//		if(random > prob)
//		{
//			fix[i] = true; // mobile
//			cont++; // counts mobile DoFs
//		}
//		else
//			fix[i] = false; // fixed
//		if(debug)
//			printf("%d ",fix[i]);
//	}
//	if(debug)
//		printf("\nmobile= %d\n",cont);
//	return cont;
//}

// Creates a random mask of fixed DiHedral coordinates according to a probability "prob"
// (Inter-chain Rotational/Translational coordinates will be always kept mobile)
// (Given a random number [0:1), an IC will be selected if "random>prob")
// "fix" array must be already allocated for full- "size" elements!
int fixRandDH(Macromolecule *mol,tri *props,bool *fix,int type,int model,float prob)
{
	bool debug = false;
	Segment *seg;
	Residue *res;
	pdbIter *iter_seg = new pdbIter( mol, true, true, true, true ); // Iterator to screen segments
	pdbIter *iter_frags; // iter residues
	int dh_index = 0;
	int ic_index = 0;
	int index_res = 0;
	int mobile = 0;
	int num_res,resn,iang_max,i,r,temp;
	int *indices=NULL; // NULL is needed by realloc!
	TMOL fragtype;

	int seg_atoms_old,interchain=0;
	int seg_atoms = 0; // should be initialized
	bool may_have_rot = false;

	if(debug)
		printf("Fixing Random (Dihedral only)\n");

	// Creating array with Dihedral coordinate indices
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = iter_seg->get_segment();
		iter_frags = new pdbIter( seg );
		num_res = iter_frags->num_fragment();

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

		// Checking single atom segments... SMOL's
		seg_atoms_old = seg_atoms;
		seg_atoms = seg->num_atoms();
		if(seg_atoms_old > 1)
			may_have_rot = true;

		// Inter-segment Internal Coordinates are allways mobile
		if(iter_seg->pos_segment !=0)
		{
			if(seg_atoms != 1 && (seg_atoms_old != 1 || may_have_rot) ) // NON single atom segments
				interchain = 6;
			else
				interchain = 3;

			if(debug)
				printf("fixRandDH> Inter-segment Internal Coords. set mobile [%d-%d]\n",ic_index,ic_index+5);
			for(i=0; i<interchain; i++,ic_index++)
			{
				fix[ic_index] = true; // mobile
				mobile++; // count mobile coordinates
			}
		}

		// Dihedral coordinates
		for( iter_frags->pos_fragment = 0; !iter_frags->gend_fragment(); iter_frags->next_fragment() )
		{
			res = ( Residue * ) iter_frags->get_fragment();
//			fragtype = res->getMolType();
			if(fragtype != tmol_smol)
				resn = resnum_from_resname( res->getName() );
			else
				resn = 666; // SMOL's can not be tabulated...

			if(debug)
				printf("fixRandDH> res= %d  fragtype= %d  resn= %d\n",iter_frags->pos_fragment,fragtype,resn);

			// if Protein fragment
			if( fragtype == tmol_protein )
			{
				// PHI
				if( !(iter_frags->pos_fragment == 0 || resn == PRO
						|| (iter_frags->pos_fragment == num_res-1 && (model == 0 || model == 3) ) ) )
				{
					// resizes indices list to allocate new DiHedral coordinate
					indices = ( int * ) realloc( indices, (dh_index+1) * sizeof( int ) );
					if(!indices)
					{
						fprintf(stderr,"Realloc failed! Forcing exit!\n");
						exit(1);
					}
					indices[dh_index] = ic_index; // store the IC index corresponding to "fix" array
					ic_index++;
					dh_index++;
				}

				// CHI
				if( type == 2 && model != 0 && model != 3)
					if( props[index_res].nan==3 ||
							(props[index_res].nan==2 &&
									(iter_frags->pos_fragment==0 ||
											(iter_frags->pos_fragment==num_res-1 && model==1) ) ) )
					{
						// resizes indices list to allocate new DiHedral coordinate
						indices = ( int * ) realloc( indices, (dh_index+1) * sizeof( int ) );
						if(!indices)
						{
							fprintf(stderr,"Realloc failed! Forcing exit!\n");
							exit(1);
						}
						indices[dh_index] = ic_index; // store the IC index corresponding to "fix" array
						ic_index++;
						dh_index++;
					}

				// PSI
				if( !(iter_frags->pos_fragment == num_res-1
						|| (iter_frags->pos_fragment == 0 && (model==0 || model==3) ) )
						|| model==2 ) // Full-Atom allways has PSI
				{
					// resizes indices list to allocate new DiHedral coordinate
					indices = ( int * ) realloc( indices, (dh_index+1) * sizeof( int ) );
					if(!indices)
					{
						fprintf(stderr,"Realloc failed! Forcing exit!\n");
						exit(1);
					}
					indices[dh_index] = ic_index; // store the IC index corresponding to "fix" array
					ic_index++;
					dh_index++;
				}
			}
			// if RNA/DNA fragment
			else if( resn == GUA || resn == ADE || resn == CYT || resn == URA ||
				resn == DGUA || resn == DADE || resn == DCYT || resn == DTHY )
			{
				if(iter_frags->pos_fragment != num_res-1)
					iang_max=6; // not last fragment
				else
					iang_max=4; // last fragment

				if(debug)
				{
					printf("fixRandDH> NAcid dihedral variables: ");
					if( type != 2 ) // if has not CHI)
						printf("%d\n",iang_max-1);
					else
						printf("%d\n",iang_max);
				}
				for(int iang=0; iang<iang_max; iang++)
				{
					if( !(iang == 3 && type != 2) ) // if not (is CHI and type!=2)
					{
						// resizes indices list to allocate new DiHedral coordinate
						indices = ( int * ) realloc( indices, (dh_index+1) * sizeof( int ) );
						if(!indices)
						{
							fprintf(stderr,"Realloc failed! Forcing exit!\n");
							exit(1);
						}
						indices[dh_index] = ic_index; // store the IC index corresponding to "fix" array
						ic_index++;
						dh_index++;
					}
				}
			}
			// For SMOL - Small MOLecules (ligands) --> Nothing "flexible" should be done here.
			else if( fragtype != tmol_smol )
			{
			    // NOT FOUND MOL-TYPE
				printf("fixRandDH> Sorry, unknown fragment type (TMOL)! Forcing exit!\n");
				exit(2);
			}

			index_res++;
		}
		delete iter_frags;
	}
	iter_seg->clean_virtual();
	delete(iter_seg);

	if(debug)
	{
		printf("fixRandDH> Dihedrals found: %d\nfixRandDH> ",dh_index);
		for(int i=0; i<dh_index; i++)
			printf("%d ",indices[i]);
		printf("\n");
	}

	// Now shuffle
	for(i=0; i<dh_index; i++)
	{
		// BUG CHUNGO QUE TE CAGAS! (17/5/2011)
		// r = dh_index * (float) rg->Random(); // [0:dh_index-1) Playing dice with Mersenne!
		r = dh_index * rg->Random(); // [0:dh_index-1) Playing dice with Mersenne!

		// swaping indices between "i" and "r"
		temp = indices[r];
		indices[r] = indices[i];
		indices[i] = temp;
	}

	// Showing shuffle
	if(debug)
	{
		fprintf(stderr,"fixRandDH> Current Shuffle: ");
		for(int i=0; i<dh_index; i++)
			fprintf(stderr,"%d ",indices[i]);
		fprintf(stderr,"\n");
	}

	// Once shuffled... the first "prob-percent" dihedral coordinates will be fixed
	for(i=0; i<dh_index; i++)
	{
		if(i < (int)((float)dh_index*prob) )
			fix[ indices[i] ] = 0; // fixed
		else
		{
			fix[ indices[i] ] = 1; // mobile
			mobile++; // count mobile coordinates
		}
	}
	free(indices);

	// Showing mask
	if(debug)
	{
		fprintf(stderr,"fixRandDH> Current Mask: ");
		for(int i=0; i<ic_index; i++)
			fprintf(stderr,"%d ",fix[i]);
		fprintf(stderr,"\n");
	}

	return mobile;
}

// CHECK THIS WITH "RNA/DNA/SMOL's"............
// Creates a mask of fixed DiHedral coordinates according to Secondary Structure
// (Inter-chain Rotational/Translational coordinates will be always kept mobile)
// "fix" array must be already allocated for full- "size" elements!
// "ss_table" is a num_res sized array with the molecule single char SS identifiers.
// "fix_ss" is an string with the residue char SS identifiers which will be fixed.
int fixSS(Macromolecule *mol,tri *props,bool *fix,char *ss_table,char *fix_ss,int type,int model)
{
	bool debug = false;
	Segment *seg;
	Residue *res;
	pdbIter *iter_seg = new pdbIter( mol ); // Iterator to screen segments
	pdbIter *iter_frags; // iter residues
	int ic_index = 0;
	int index_res = 0;
	int mobile = 0;
	int num_res,resn,iang_max,i;
	TMOL fragtype;
	int ss_len = strlen(fix_ss); // gets number of SS identifiers
	bool fixit;

	// Screening Macromolecule
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = iter_seg->get_segment();
		iter_frags = new pdbIter( seg );
		num_res = iter_frags->num_fragment();

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

		// Inter-segment Internal Coordinates are allways mobile
		if(iter_seg->pos_segment !=0)
		{
			if(debug)
				printf("fixSS> Inter-segment Internal Coords. set mobile [%d-%d]\n",ic_index,ic_index+5);
			for(i=0; i<6; i++,ic_index++)
			{
				fix[ic_index] = true; // mobile
				mobile++; // count mobile coordinates
			}
		}

		// Dihedral coordinates
		for( iter_frags->pos_fragment = 0; !iter_frags->gend_fragment(); iter_frags->next_fragment() )
		{
			res = ( Residue * ) iter_frags->get_fragment();
			resn = resnum_from_resname( res->getName() );
//			fragtype = res->getMolType();

			fixit = false; // by default mobile
			for( i=0; i < ss_len; i++)
				if(fix_ss[i] == ss_table[index_res]) // if any SS id matches...
					fixit = true;

			// if Protein fragment
			if( fragtype == tmol_protein )
			{
				// PHI
				if( !(iter_frags->pos_fragment == 0 || resn == PRO
						|| (iter_frags->pos_fragment == num_res-1 && (model == 0 || model == 3) ) ) )
				{
					if(fixit)
						fix[ic_index] = false; // fixed
					else
					{
						fix[ic_index] = true; // mobile
						mobile++; // count mobile coordinates
					}
					ic_index++;
				}

				// CHI
				if( type == 2 && model != 0 && model != 3)
					if( props[index_res].nan==3 ||
							(props[index_res].nan==2 &&
									(iter_frags->pos_fragment==0 ||
											(iter_frags->pos_fragment==num_res-1 && model==1) ) ) )
					{
						if(fixit)
							fix[ic_index] = false; // fixed
						else
						{
							fix[ic_index] = true; // mobile
							mobile++; // count mobile coordinates
						}
						ic_index++;
					}

				// PSI
				if( !(iter_frags->pos_fragment == num_res-1
						|| (iter_frags->pos_fragment == 0 && (model==0 || model==3) ) )
						|| model==2 ) // Full-Atom allways has PSI
				{
					if(fixit)
						fix[ic_index] = false; // fixed
					else
					{
						fix[ic_index] = true; // mobile
						mobile++; // count mobile coordinates
					}
					ic_index++;
				}
			}
			// if RNA fragment
			else if( resn == GUA || resn == ADE || resn == CYT || resn == URA )
			{
				if(iter_frags->pos_fragment != num_res-1)
					iang_max=6; // not last fragment
				else
					iang_max=4; // last fragment

				if(debug)
				{
					printf("fixSS> NAcid dihedral variables: ");
					if( type != 2 ) // if has not CHI)
						printf("%d\n",iang_max-1);
					else
						printf("%d\n",iang_max);
				}
				for(int iang=0; iang<iang_max; iang++)
				{
					if( !(iang == 3 && type != 2) ) // if not (is CHI and type!=2)
					{
						if(fixit)
							fix[ic_index] = false; // fixed
						else
						{
							fix[ic_index] = true; // mobile
							mobile++; // count mobile coordinates
						}
						ic_index++;
					}
				}
			}
			else
			{
				printf("fixSS> Sorry, unknown fragment type! Forcing exit!\n");
				exit(1);
			}

			index_res++;
		}
		delete iter_frags;
	}
	delete(iter_seg);

	// Showing mask
	if(debug)
	{
		printf("fixSS> Current Mask: ");
		for(int i=0; i<ic_index; i++)
			printf("%d ",fix[i]);
		printf("\nmobile= %d\n",mobile);
	}
	return mobile;
}

// Extracts Dihedral-Angle components directly from the Hessian Matrix raw modes
// Allocates the mode array, only for the first time.
// (Please, free it when not needed)
//void dihedral_comps(Macromolecule *mol,double *evec,int nev,int size,tri *props,double **pp_mode)
void dihedral_comps(double *evec,int nev,int size,double **pp_mode)
{
	int offset;
	pdbIter *iter;
	Residue *res;
	double *p_mode = *pp_mode;

	if(p_mode == NULL) // then, allocate mem.
	{
		// Allocates more elements than needed!
		if(	!(p_mode = (double *) malloc( sizeof(double) * size ) ) )
		{
			printf("Msg(dihedral_comps): mode array not allocated!\nForcing Exit!\n");
			exit(1);
		}
		*pp_mode = p_mode; // outputs mode allocated address
	}

	offset = nev * size;
	for(int k=0; k<size; k++)
		p_mode[k] = evec[ offset + k ]; // retrieving mode
}

// Change normal modes (internal coordinates) from a "i" model into a "f" model.
// WARNINGS: 1) Proteins only! (no CG-models in DNA/RNA/SMOL)
//           2) "i" is ALLWAYS the finer-grained model, "f" is ALLWAYS the coarser-grained model.
// inverse = false (default): Matching "f" DoFs will be set to "i"'s ones (all "f" DoFs will be filled)
// inverse = true: Matching "i" DoFs will be set to "f"'s ones, (not-matching "i" DoFs will be = zero)
void changemodelIC(Macromolecule *moli, tri *propsi, tri *propsf, floating *eveci, floating *evecf, bool *fixi, bool *fixf, int nmodes, int sizei, int sizef, int model, int type, int modelf, int typef, bool inverse)
{
	bool debug = false;
	int j = 0; // "i" current dihedral index (screens all dihedrals)
	int jjf = 0; // "f" current dihedral index (screens all dihedrals)
	int ji = 0; // "i" active dihedral index (screens all mobile dihedrals)
	int jf = 0; // "f" active dihedral index (screens all mobile dihedrals)
	int num_res,resn,res_index;

	pdbIter *iter_seg,*iter_frag;
	iter_seg = new pdbIter( moli ); // Iterator to screen segments "i"

	Segment *seg;
	Residue *res;
	Atom *atom;
	TMOL fragtype;

	// Screening segments
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

		if(debug)
			printf("\nProcessing segment %d (%d): \n",iter_seg->pos_segment,num_res);

		if(iter_seg->pos_segment != 0) // non-first segment
		{
			// Check whether any of the 6D inter-segment variables are fixed
			// 3 TRANSLATIONS + 3 ROTATIONS
			for(int axis=0; axis<6; axis++)
			{
				if( fixi == NULL || fixi[j] ) // Check whether current variable it's fixed
				{
					if(inverse) // inverse
						for(int n=0; n<nmodes; n++) // screening only modes to save
							eveci[n*sizei+ji] = evecf[n*sizef+jf];
					else // normal
						for(int n=0; n<nmodes; n++) // screening only modes to save
							evecf[n*sizef+jf] = eveci[n*sizei+ji];
					ji++;
					jf++; // both will have the same
				}
				j++; // "i" absolute index
				jjf++; // "f" absolute index
			}
		}

		// Screen ALL fragments
		for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			res = ( Residue * ) iter_frag->get_fragment();
			resn = resnum_from_resname( res->getName() );
//			fragtype = res->getMolType();

			if(debug)
				printf (" Processing fragment %d %s (%d) type= %d\n",iter_frag->pos_fragment,res->getName(),resn,fragtype);

			// if Protein fragment
			if( fragtype == tmol_protein )
			{
				// ********************************************************************************************
				// PHI --> NOT FIRST, NOT PRO
				if ( ( model != 0 && !(iter_frag->pos_fragment == 0 || resn == PRO ) ) ||
				     ( model == 0 && !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) ) )
				{
					//  printf (" Not first not PRO -> PHI \n");
					if( fixi == NULL || fixi[ j ] ) // Check whether current variable it's fixed
					{
						if(modelf == 0) // CA-model
						{
							// PHI --> NOT FIRST, NOT PRO, NOT LAST too! (CA-model)
							if( !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) )
							{ // "f" has got PHI
								if(inverse) // inverse
									for(int n=0; n<nmodes; n++) // screening only modes to save
										eveci[n*sizei+ji] = evecf[n*sizef+jf];
								else // normal
									for(int n=0; n<nmodes; n++) // screening only modes to save
										evecf[n*sizef+jf] = eveci[n*sizei+ji];
								jf++;
								jjf++; // "f" absolute index
							}
							else // "f" hasn't got PHI
								if(inverse)
									for(int n=0; n<nmodes; n++) // screening only modes to save
										eveci[n*sizei+ji] = 0.0; // setting zeros into "i"
						}
						else // 3BB2R and Full-Atom share PHIs
						{
							if(inverse) // inverse
								for(int n=0; n<nmodes; n++) // screening only modes to save
									eveci[n*sizei+ji] = evecf[n*sizef+jf];
							else // normal
								for(int n=0; n<nmodes; n++) // screening only modes to save
									evecf[n*sizef+jf] = eveci[n*sizei+ji];
							jf++;
							jjf++; // "f" absolute index
						}
						ji++;
					}
					else
					{
						if(modelf == 0) // CA-model
						{
							// PHI --> NOT FIRST, NOT PRO, NOT LAST too! (CA-model)
							if( !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) )
							{ // "f" has got PHI
								jjf++; // "f" absolute index
							}
						}
						else // 3BB2R and Full-Atom share PHIs
						{
							jjf++; // "f" absolute index
						}
					}
					j++; // "i" index
				}  // NOT FIRST NOT PRO

				// ********************************************************************************************
				//  LATERAL CHAIN-->CHI
				// 3 dihedrals (normal residue) or 2 dihedrals (ending residue)
				if(type == 2) // phi,chi,psi,...
					if( propsi[res_index].nan==3 ||
							(propsi[res_index].nan==2 &&
											(iter_frag->pos_fragment==0 ||
													(iter_frag->pos_fragment==num_res-1 && model==1) ) ) )
					{
						//printf (" Lateral Chain --> CHI\n");
						if( fixi == NULL || fixi[ j ] ) // Check whether the variable it's fixed
						{
							if(typef == 2 && modelf != 0 ) // if "f" has CHIs
							{
								if(inverse) // inverse
									for(int n=0; n<nmodes; n++) // screening only modes to save
										eveci[n*sizei+ji] = evecf[n*sizef+jf];
								else // normal
					 				for(int n=0; n<nmodes; n++) // screening only modes to save
										evecf[n*sizef+jf] = eveci[n*sizei+ji];
								jf++; // both will have the same
								jjf++; // "f" absolute index
							}
							else // "f" hasn't got CHI
								if(inverse)
									for(int n=0; n<nmodes; n++) // screening only modes to save
										eveci[n*sizei+ji] = 0.0; // setting zeros into "i"
							ji++;
						}
						j++; // der[] index
					}

				// Just for 20 --> 12 (normal) and 12 --> 20 (inverse) cases
				if(model == 2 && type == 0 && modelf == 1 && typef == 2) // 20 --> 12 (normal) and 12 --> 20 (inverse) cases
					if( propsf[res_index].nan==3 ||
							(propsf[res_index].nan==2 &&
											(iter_frag->pos_fragment==0 ||
													(iter_frag->pos_fragment==num_res-1 && modelf==1) ) ) )
					{
						//printf (" Lateral Chain --> CHI\n");
						if( fixf == NULL || fixf[ jjf ] ) // Check whether the variable it's fixed
						{
							// "f" has CHIs (see above)
							if(!inverse) // normal
								for(int n=0; n<nmodes; n++) // screening only modes to save
									evecf[n*sizef+jf] = 0.0; // setting zeros into "f"
							jf++; // only "f" has got CHI
						}
						jjf++; // "f" has got CHI ("i" hasn't got it... so no j++)
					}

				// ********************************************************************************************
				// NOT LAST RESIDUE--> PSI
				//
				// "PSI-bond" (in Full-Atom, every residue has PSI)
				if( ( iter_frag->pos_fragment != num_res - 1 && model == 1 ) || model==2 || // Full-Atom allways has PSI
				    ( iter_frag->pos_fragment != 0 && iter_frag->pos_fragment != num_res-1 && model == 0) ) // CA-case
				{
					if( fixi == NULL || fixi[ j ] ) // Check whether the variable it's fixed
					{
						if( modelf == 0) // CA-model
						{
							// "PSI-bond" (In CA-only model, non-first and non-last)
							if( !(iter_frag->pos_fragment == num_res-1 || iter_frag->pos_fragment == 0) ) // if "f" has got PSI
							{
								if(inverse) // inverse
									for(int n=0; n<nmodes; n++) // screening only modes to save
										eveci[n*sizei+ji] = evecf[n*sizef+jf];
								else // normal
									for(int n=0; n<nmodes; n++) // screening only modes to save
										evecf[n*sizef+jf] = eveci[n*sizei+ji];
								jf++; // both will have the same
								jjf++; // "f" absolute index
							}
							else // "f" hasn't got PSI
								if(inverse)
									for(int n=0; n<nmodes; n++) // screening only modes to save
										eveci[n*sizei+ji] = 0.0; // setting zeros into "i"
						}
						else // 3BB2R and Full-Atom
							if ( iter_frag->pos_fragment != num_res - 1 || modelf==2 ) // if "f" has got PSI
							{
								if(inverse) // inverse
									for(int n=0; n<nmodes; n++) // screening only modes to save
										eveci[n*sizei+ji] = evecf[n*sizef+jf];
								else // normal
									for(int n=0; n<nmodes; n++) // screening only modes to save
										evecf[n*sizef+jf] = eveci[n*sizei+ji];
								jf++; // both will have the same
								jjf++; // "f" absolute index
							}
							else // "f" hasn't got PSI
								if(inverse)
									for(int n=0; n<nmodes; n++) // screening only modes to save
										eveci[n*sizei+ji] = 0.0; // setting zeros into "i"
						ji++;
					}
					else
					{
						if(modelf == 0) // CA-model
						{
							// "PSI-bond" (In CA-only model, non-first and non-last)
							if( !(iter_frag->pos_fragment == num_res-1 || iter_frag->pos_fragment == 0) ) // if "f" has got PSI
							{ // "f" has got PSI
								jjf++; // "f" absolute index
							}
						}
						else // 3BB2R and Full-Atom share PHIs
						{
							if ( iter_frag->pos_fragment != num_res - 1 || modelf==2 ) // if "f" has got PSI
								jjf++; // "f" absolute index
						}
					}
					j++; // der[] index
				}
			}
			else // NOT FOUND MOL-TYPE
			{
				printf("Msg(changemodelIC): Given CG-models are only available for proteins,\n"
						"different CG models are senseless in mixed structures with DNA/RNA/SMOL\n"
						"Please, use only protein structures with this option. Forcing exit!\n");
				exit(2);
			}

			res_index++;
		}
	} // end derivatives
}

// Change fix-arrays (internal coordinates) from a "i" model into a "f" model.
// WARNINGS: 1) Proteins only! (no CG-models in DNA/RNA/SMOL)
//           2) "i" is ALLWAYS the finer-grained model, "f" is ALLWAYS the coarser-grained model.
//             (no memory allocation is performed)
// inverse = false (default): Matching "f" DoFs will be set to "i"'s ones (all "f" DoFs will be filled)
// inverse = true: Matching "i" DoFs will be set to "f"'s ones, (not-matching "i" DoFs will be = zero)
void changefixIC(Macromolecule *moli, tri *propsi, tri *propsf, bool *fixi, bool *fixf,
		int model, int type, int modelf, int typef, bool inverse)
{
	bool debug = false;
	int j = 0; // "i" current dihedral index (screens all dihedrals)
	int jjf = 0; // "f" current dihedral index (screens all dihedrals)
	int num_res,resn,res_index;

	if(inverse) // inverse
	{
		if(fixf == NULL) // nothing to do
		{
			fixi = NULL;
			return; // exit
		}
	}
	else // normal
		if(fixi == NULL) // nothing to do
		{
			fixf = NULL;
			return; // exit
		}

	pdbIter *iter_seg,*iter_frag;
	iter_seg = new pdbIter( moli ); // Iterator to screen segments "i"

	Segment *seg;
	Residue *res;
	Atom *atom;
	TMOL fragtype;

	// Screening segments
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

		if(debug)
			printf("\nProcessing segment %d (%d): \n",iter_seg->pos_segment,num_res);

		if(iter_seg->pos_segment != 0) // non-first segment
			// Check whether any of the 6D inter-segment variables are fixed
			// 3 TRANSLATIONS + 3 ROTATIONS  // both should have the same, allways
			for(int axis=0; axis<6; axis++)
			{
				if(inverse)
					fixi[j] = fixf[jjf];
				else
					fixf[jjf] = fixi[j];
				j++; // "i" index
				jjf++; // "f" index
			}

		// Screen ALL fragments
		for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			res = ( Residue * ) iter_frag->get_fragment();
			resn = resnum_from_resname( res->getName() );
//			fragtype = res->getMolType();

			if(debug)
				printf (" Processing fragment %d %s (%d) type= %d\n",iter_frag->pos_fragment,res->getName(),resn,fragtype);

			// if Protein fragment
			if( fragtype == tmol_protein )
			{
				// ********************************************************************************************
				// PHI --> NOT FIRST, NOT PRO
				// "i" allways has got PHI (model==1,2), excepting model==0
				if ( ( model != 0 && !(iter_frag->pos_fragment == 0 || resn == PRO ) ) ||
				     ( model == 0 && !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) ) )
				{
					//  printf (" Not first not PRO -> PHI \n");
					if(modelf == 0) // "f" CA-model
					{
						// PHI --> NOT FIRST, NOT PRO, NOT LAST too! (CA-model)
						if( !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) )
						{ // "f" has got PHI
							if(inverse)
								fixi[j] = fixf[jjf];
							else
								fixf[jjf] = fixi[j];
							jjf++; // "f" index
						}
						else // "f" hasn't got PHI
							if(inverse)
								fixi[j] = false; // "i"'s PHI is fixed
					}
					else // 3BB2R and Full-Atom share PHIs
					{
						if(inverse)
							fixi[j] = fixf[jjf];
						else
							fixf[jjf] = fixi[j];
						jjf++; // "f" index
					}
					j++; // "i" index
				}  // NOT FIRST NOT PRO

				// ********************************************************************************************
				//  LATERAL CHAIN-->CHI
				// 3 dihedrals (normal residue) or 2 dihedrals (ending residue)
				if(type == 2) // phi,chi,psi,...
					if( propsi[res_index].nan==3 ||
							(propsi[res_index].nan==2 &&
											(iter_frag->pos_fragment==0 ||
													(iter_frag->pos_fragment==num_res-1 && model==1) ) ) )
					{
						if(typef == 2 && modelf != 0 ) // if "f" has CHIs
						{
							if(inverse)
								fixi[j] = fixf[jjf];
							else
								fixf[jjf] = fixi[j];
							jjf++; // "f" index
						}
						else // "f" hasn't got CHI
							if(inverse)
								fixi[j] = false; // "i"'s CHI is fixed

						j++; // "i" index ("i" has got CHI)
					}

				// if "f" has got CHI but "i" hasn't.
				// Just for 20 --> 12 (normal) and 12 --> 20 (inverse) cases
				if(model == 2 && type == 0 && modelf == 1 && typef == 2) // 20 --> 12 (normal) and 12 --> 20 (inverse) cases
					if( propsf[res_index].nan==3 ||
							(propsf[res_index].nan==2 &&
											(iter_frag->pos_fragment==0 ||
													(iter_frag->pos_fragment==num_res-1 && modelf==1) ) ) )
					{
						//printf (" Lateral Chain --> CHI\n");
						// "f" has CHIs (see above)
						if(!inverse) // normal
							fixf[jjf] = false; // "f"-CHIs fixed because "i" hasn't got it.
						jjf++; // "f" has got CHI ("i" hasn't got it... so no j++)
					}

				// ********************************************************************************************
				// NOT LAST RESIDUE--> PSI
				// "PSI-bond" (in Full-Atom, every residue has PSI)
				if( ( iter_frag->pos_fragment != num_res - 1 && model == 1 ) || model==2 || // Full-Atom allways has PSI
				    ( iter_frag->pos_fragment != 0 && iter_frag->pos_fragment != num_res-1 && model == 0) ) // CA-case
				{ // if "i" has got PSI
					if( modelf == 0 ) // CA-model
					{
						// "PSI-bond" (In CA-only model, non-first and non-last)
						if( !(iter_frag->pos_fragment == num_res-1 || iter_frag->pos_fragment == 0) )
						{ // if "f" has got PSI
							if(inverse)
								fixi[j] = fixf[jjf];
							else
								fixf[jjf] = fixi[j];
							jjf++; // "f" index
						}
						else // "f" hasn't got PSI
							if(inverse)
								fixi[j] = false; // "i"'s PHI is fixed
					}
					else // 3BB2R and Full-Atom models
					{
						if ( iter_frag->pos_fragment != num_res - 1 || modelf==2 )
						{ // if "f" has got PSI
							if(inverse)
								fixi[j] = fixf[jjf];
							else
								fixf[jjf] = fixi[j];
							jjf++; // "f" index
						}
						else // "f" hasn't got PSI
							if(inverse)
								fixi[j] = false; // "i"'s PHI is fixed
					}
					j++; // der[] index
				}
			}
			else // NOT FOUND MOL-TYPE
			{
				printf("Msg(changefixIC): Given CG-models are only available for proteins,\n"
						"different CG models are senseless in mixed structures with DNA/RNA/SMOL\n"
						"Please, use only protein structures with this option. Forcing exit!\n");
				exit(2);
			}
			res_index++;
		}
	} // end derivatives

	if(debug)
		printf("Msg(changefixIC):  i --> %d DoFs   f --> %d DoFs\n",j,jjf);
}

// Pablo's "inverse exponential function"
double inv_exp(double k, double x, double x0, double power)
{
	return ( k / ( 1.0 + pow( x/x0, power ) ) );
}

// Normalices eigenvectors
int norm_evec(double *evec,int nevec,int size)
{
	double norm = 0.0;
	int offset;
	// Normalizing eigenvectors
	for(int i=0; i<nevec; i++)
	{
		norm = 0.0;
		offset = i * size;
		for(int j=0; j<size; j++) norm += pow(evec[offset+j],2); // computes the norm
		norm = sqrt(norm);
		for(int j=0; j<size; j++) evec[offset+j] /= norm ; // this normalizes
	}
	return( 0 ); // normal exit status
}

// Selects the "n_chain"th Chain from a Macromolecule
// *****   Watch out with NUCLEIC ACIDS and SMOL's ******
Macromolecule *select_chain(Macromolecule *in, int n_chain)
{
	Macromolecule *out;
	Protein *pr2;
	Chain *ch;
	pdbIter *it_ch;
	out = new Macromolecule( "selected chain" );
	pr2 = new Protein();
	it_ch = new pdbIter( in );
	it_ch->pos_chain = n_chain; // "n_chain" chain index
	ch = it_ch->get_chain(); // get first chain
	pr2->add( ch ); // adding chain
	out->add( pr2 ); // adding protein
	delete it_ch;
	return out;
}

// Selects the "n_seg"th Segment from a Macromolecule
// *****   Watch out with NUCLEIC ACIDS and SMOL's ******
Macromolecule *select_segment(Macromolecule *in, int n_seg)
{
	Macromolecule *out;
	Protein *pr2;
	Chain *ch2;
	Segment *seg;
	pdbIter *it_seg;
	out = new Macromolecule( "selected chain" );
	pr2 = new Protein();
	ch2 = new Chain();
	it_seg = new pdbIter( in, true, true, true, true );
	it_seg->pos_segment = n_seg; // "n_seg" segment index
	seg = it_seg->get_segment(); // get segment
	ch2->add( seg ); // adding segment
	pr2->add( ch2 ); // adding chain
	out->add( pr2 ); // adding protein
	delete it_seg;
	return out;
}

//// Selects the "n_seg"th Segment from a Macromolecule
//// *****   Watch out with NUCLEIC ACIDS and SMOL's ******
//Macromolecule *select_segment(Macromolecule *in, int n_seg)
//{
//	Macromolecule *out;
//	Protein *pr2;
//	Chain *ch2;
//	Segment *seg;
//	pdbIter *it_seg;
//	out = new Macromolecule( "selected chain" );
//	pr2 = new Protein();
//	ch2 = new Chain();
//	it_seg = new pdbIter( in );
//	it_seg->pos_segment = n_seg; // "n_seg" segment index
//	seg = it_seg->get_segment(); // get segment
//	ch2->add( seg ); // adding segment
//	pr2->add( ch2 ); // adding chain
//	out->add( pr2 ); // adding protein
//	delete it_seg;
//	return out;
//}

// Number of Degrees of Freedom in a Macromolecule "mol" (given "props")
int num_dofs(Macromolecule *mol,tri *props, int *size_dh, int *n_chain, int *n_seg)
{
	bool verb=false;
	pdbIter *iter;
	int size=0;
	int seg=0,chain;
	int old_size;

	// Computing Total number of degrees of freedom (hessian matrix rank)
	iter = new pdbIter( mol );
	for ( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
		size += props[iter->pos_fragment].nan; // Each residue may have different number of dihedral-angles!

	old_size = size; // temp
	seg = iter->num_segment();
	size += 6*(seg-1); // (n_seg-1)*6 aditional degrees of freedom must be added
	chain = iter->num_chain();
	delete iter;

	if(verb)
	{
		printf( "Msg(num_dofs): Number of dihedral angles: %d\n", old_size );
		printf( "Msg(num_dofs): Number of Inter-segment coords: %d (Rot+Trans)\n", size-old_size );
		printf( "Msg(num_dofs): Number of Internal Coordinates: %d (Hessian rank)\n", size );
		printf( "Msg(num_dofs): Number of Chains: %d   Segments: %d\n", chain, seg );
	}

	if(size_dh!=NULL)
		*size_dh = old_size;
	if(n_chain!=NULL)
		*n_chain = chain;
	if(n_seg!=NULL)
		*n_seg = seg;

	return size;
}

//// Checking whether memory was properly allocated in "ptr" by malloc, and showing "info"
//void check_pointer(void *ptr,char *info)
//{
//	if(!ptr)
//	{
//		fprintf(stderr,"Memory allocation failed!\n");
//		fprintf(stderr,info);
//		fprintf(stderr,"\nForcing exit from check_pointer()\n");
//		exit(1);
//	}
//}

// Fills a vector with random values from "low" to "high" [low:high)
void rand_vector(double *v, int size, double low, double high)
{
	int i;
	double range;
	range = high - low;

	for(i=0; i<size; i++)
		v[i] = low + range * (double) rg->Random(); // [0:1) Playing dice with Mersenne!
}

// Sort the "ipas" to feed directly the Hessian matrix computation routines.
// (Quick-sort implementation for "ipas")
twid **sort_ipas(twid *ipas, int *pnipas, int *unat, int num_units)
{
 twid **pipas; // Pointer to "ipas" array
 int i,x;
 int nremoved = 0;
 int nipas;

 nipas = *pnipas;  // Output new number of nipas (due to removed intra-unit interactions)

 // Allocating a "pointer to ipas" array (this array will be sorted later...)
 if( !(pipas = (twid **) malloc( sizeof(twid *) * nipas ) ) ) // Array of pointers to "ipas"
 {
	 printf("Msg(sort_ipas): Sorry, memory allocation failed!\n");
	 exit(1);
 }

 // Loading "pointer to ipas" array
// printf("Loading ipas (nipas= %d):\n", nipas);
 for(i = 0; i < nipas; i++)
 {
	 x = i-nremoved;
	 pipas[x] = ipas + i; // pointer arithmetics
	 pipas[x]->i = unat[pipas[x]->k]; // Filling units indices (i,j)
	 pipas[x]->j = unat[pipas[x]->l]; // Filling units indices (i,j)

//	 fprintf(stderr,"i= %d  nipas= %d  x= %d  pipas[x]->i= %d  pipas[x]->j= %d  unat[pipas[x]->k]= %d unat[pipas[x]->l]= %d\n",i,nipas,x,pipas[x]->i,pipas[x]->j,unat[pipas[x]->k],unat[pipas[x]->l]);

//	 // Easy way to get num_units
//	 if(pipas[x]->i > num_units)
//		 num_units = pipas[x]->i;
//	 if(pipas[x]->j > num_units)
//		 num_units = pipas[x]->j;

	 if(pipas[x]->i == pipas[x]->j)
		 nremoved++; // This removes Intra-unit contacts!

	 // printf("ipa= %d  Un(k,l)= %d,%d  k,l= %d,%d  i,j= %d,%d\n",i,unat[pipas[x]->k],unat[pipas[x]->l],pipas[x]->k,pipas[x]->l,pipas[x]->i,pipas[x]->j);

//	 if(pipas[x]->i != pipas[x]->j)
//		 printf("ipa= %d  Un(k,l)= %d,%d  k,l= %d,%d  i,j= %d,%d\n",i,unat[pipas[x]->k],unat[pipas[x]->l],pipas[x]->k,pipas[x]->l,pipas[x]->i,pipas[x]->j);
//	 else
//	 {
//		 printf("ipa= %d  Un(k,l)= %d,%d  k,l= %d,%d  i,j= %d,%d REMOVED!\n",i,unat[pipas[x]->k],unat[pipas[x]->l],pipas[x]->k,pipas[x]->l,pipas[x]->i,pipas[x]->j);
//		 nremoved++; // This removes Intra-unit contacts!
//	 }

 }
 nipas -= nremoved; // Intra-unit contacts removed!

// // Print the original array
// printf("Before quicksort:\n");
// for(i = 0; i < nipas; i++)
//	 printf("ipa= %d  unat(k,l)= %d,%d  k,l= %d,%d  i,j= %d,%d\n",i,unat[pipas[i]->k],unat[pipas[i]->l],pipas[i]->k,pipas[i]->l,pipas[i]->i,pipas[i]->j);
// printf("\n");

 // quicksort(array, 0, (MAXARRAY - 1));
// quicksort(pipas, (long int)0,(long int) (nipas-1), num_units);
 quicksort(pipas, 0, nipas-1, num_units);

// // print the "quicksorted" array
// fprintf(stderr,"After  quicksort:\n");
// for(i = 0; i < nipas; i++)
//	 fprintf(stderr,"ipa= %d  Un(k,l)= %d,%d  k,l= %d,%d  i,j= %d,%d\n",i,unat[pipas[i]->k],unat[pipas[i]->l],pipas[i]->k,pipas[i]->l,pipas[i]->i,pipas[i]->j);
// fprintf(stderr,"\n");
// fprintf(stderr,"\nEND sort_ipas()\n");

 *pnipas = nipas; // Output new number of nipas (due to removed intra-unit interactions)
 return(pipas);
// exit(0);
}

/* sort everything in between `low' <-> `high' */
// Quick-sort recursive routine
void quicksort0(twid **arr, int low, int high, int size)
{
	int i = low;
	int j = high;
//	int y = 0;

	/* compare value */
//	int z = arr[(low + high) / 2];
	twid *z = arr[(low + high) / 2]; // i + N*(N-1)/2 - j*(j+1)/2;
	twid *y;
//	int max = size*(size-1)/2; // off diagonal
	int max = size*(size+1)/2; // off+on diagonal
//	int z_point = z->i + max - z->j*(z->j+1)/2;
	int z_point = z->i + max - (z->j+1)*(z->j+2)/2;

//	int indice;
//	int N=10;
//	printf("Testing index:\n");
//	for(int i=0; i<10; i++) // row --> i=[0,N-1]
//	{
//		for(int j=0; j<10; j++) // col --> j=[0,N-1]
//		{
//			if(j>=i)
//			{
////				indice = i + N*(N-1)/2 - j*(j+1)/2; // off-diagonal
//				indice = i + N*(N+1)/2 - (j+1)*(j+2)/2; // off+on diagonal
//				printf("%3d ",indice);
//			}
//			else
//				printf("%3d ",0);
//		}
//		printf("\n");
//	}
//	printf("END test\n");
//	exit(0);

//	printf("quicksort: low= %d  high= %d  z_point= %d\n",low,high,z_point);
	/* partition */
	do {
		/* find member above ... */
//		while(arr[i] < z)
//		while(arr[i]->i + max - arr[i]->j*(arr[i]->j+1)/2 < z_point ) // i + N*(N-1)/2 - j*(j+1)/2;
		while(arr[i]->i + max - (arr[i]->j+1)*(arr[i]->j+2)/2 < z_point ) // i + N*(N+1)/2 - (j+1)*(j+2)/2;
			i++;

		/* find element below ... */
//		while(arr[j] > z)
//		while(arr[j]->i + max - arr[j]->j*(arr[j]->j+1)/2 > z_point ) // i + N*(N-1)/2 - j*(j+1)/2;
		while(arr[j]->i + max - (arr[j]->j+1)*(arr[j]->j+2)/2 > z_point ) // i + N*(N+1)/2 - (j+1)*(j+2)/2;
			j--;

		// printf("\ti= %d  j= %d  arri= %d  arrj= %d  ",i,j,arr[i]->i + max - arr[i]->j*(arr[i]->j+1)/2,arr[j]->i + max - arr[j]->j*(arr[j]->j+1)/2);

		if(i <= j)
		{
			// printf("swaping i by j");
			// Swap two elements
			y = arr[i];
			arr[i] = arr[j];
			arr[j] = y;
			i++;
			j--;
		}
//		else
//			printf("no-swaping");
//		printf("\n");
	} while(i <= j);

	/* recurse */
	if(low < j)
		quicksort(arr, low, j, size);

	if(i < high)
		quicksort(arr, i, high, size);
}

/* sort everything in between `low' <-> `high' */
// Quick-sort recursive routine
void quicksortGOOD(twid **arr, long low, long high, long size)
{
	long i = low;
	long j = high;
	long laii,laij,laji,lajj;
//	int y = 0;

	/* compare value */
//	int z = arr[(low + high) / 2];
//	fprintf(stderr, "(low+high)/2 = %ld\n",(low + high) / 2);
	twid *z = arr[(low + high) / 2]; // i + N*(N-1)/2 - j*(j+1)/2;
	twid *y;
//	int max = size*(size-1)/2; // off diagonal
	long max = size*(size+1)/2; // off+on diagonal
//	int z_point = z->i + max - z->j*(z->j+1)/2;
	long lzi = (long)(z->i);
	long lzj = (long)(z->j);
	long z_point = lzi + max - (lzj+1)*(lzj+2)/2;

	/* partition */
	do {
		/* find member above ... */
		laii = (long)(arr[i]->i);
		laij = (long)(arr[i]->j);
//		while(arr[i]->i + max - (arr[i]->j+1)*(arr[i]->j+2)/2 < z_point )
		while(laii + max - (laij+1)*(laij+2)/2 < z_point )
		{
			i++;
			laii = (long)(arr[i]->i);
			laij = (long)(arr[i]->j);
		}

		/* find element below ... */
		laji = (long)(arr[j]->i);
		lajj = (long)(arr[j]->j);
//		while(arr[j]->i + max - (arr[j]->j+1)*(arr[j]->j+2)/2 > z_point )
		while(laji + max - (lajj+1)*(lajj+2)/2 > z_point )
		{
			j--;
			laji = (long)(arr[j]->i);
			lajj = (long)(arr[j]->j);
		}

//		fprintf(stderr,"\ti= %d  j= %d  arri= %d  arrj= %d  ",i,j,arr[i]->i + max - arr[i]->j*(arr[i]->j+1)/2,arr[j]->i + max - arr[j]->j*(arr[j]->j+1)/2);

		if(i <= j)
		{
			// printf("swaping i by j");
			// Swap two elements
			y = arr[i];
			arr[i] = arr[j];
			arr[j] = y;
			i++;
			j--;
		}
//		else
//			fprintf(stderr,"no-swaping");
//		fprintf(stderr,"\n");
	} while(i <= j);

	/* recurse */
	if(low < j)
		quicksort(arr, low, j, size);

	if(i < high)
		quicksort(arr, i, high, size);
}


// Sorts everything in between `low' <-> `high'
// Quick-sort recursive routine (with "long int" in "max" for big systems)
void quicksort(twid **arr, int low, int high, int size)
{
	int i = low;
	int j = high;

	/* compare value */
	twid *z = arr[(low + high) / 2];
	twid *y;
	long max = (long)size*((long)size+1)/2; // off+on diagonal
//	int z_point = z->i + max - z->j*(z->j+1)/2;
	long z_point = (long)(z->i) + max - ((long)(z->j+1))*((long)(z->j+2))/2;

	/* partition */
	do {
		/* find member above ... */
//		while(arr[i] < z)
//		while(arr[i]->i + max - arr[i]->j*(arr[i]->j+1)/2 < z_point ) // i + N*(N-1)/2 - j*(j+1)/2;
		while( max + (long)(arr[i]->i) - ((long)(arr[i]->j+1))*((long)(arr[i]->j+2))/2 < z_point ) // i + N*(N+1)/2 - (j+1)*(j+2)/2;
			i++;

		/* find element below ... */
//		while(arr[j] > z)
//		while(arr[j]->i + max - arr[j]->j*(arr[j]->j+1)/2 > z_point ) // i + N*(N-1)/2 - j*(j+1)/2;
		while( max + (long)(arr[j]->i) - ((long)(arr[j]->j+1))*((long)(arr[j]->j+2))/2 > z_point ) // i + N*(N+1)/2 - (j+1)*(j+2)/2;
			j--;

		if(i <= j)
		{
			// Swap two elements
			// printf("swaping i by j");
			y = arr[i];
			arr[i] = arr[j];
			arr[j] = y;
			i++;
			j--;
		}
//		else
//			printf("no-swaping");
//		printf("\n");
	} while(i <= j);

	/* recurse */
	if(low < j)
		quicksort(arr, low, j, size);

	if(i < high)
		quicksort(arr, i, high, size);
}

//
void swapSPmatrix(char *file, long int mem)
{
// HESSIAN COMPUTATION -->	hess_matrix[a + b*(b+1)/2] = temp;
// (a=rows, b=cols)
// 	for(long int b=size-1; b >= 0; b--) // b --> dihedrals (column)
//		for(long int a=0; a <= b; a++)	 // a --> dihedrals (row)

// BINARY SP MATRIX	-->	fwrite(&matrix[i + j*(j+1)/2],sizeof(floating),1,f_file);
// (i=rows, j=cols)

	bool debug = false;
	long int block = mem;
	long int i,j,a,b,x,y,m,c,d,next;
	long int n=0; // n --> number of elements given "block" and "size"
	long int size;
	long int stripe=0;
	long int raiz=0;
	int p=0; // p --> number of rows given "block" and "size"
	int r=0; // r --> number of remaining rows
	FILE *f_file,*f_out;
	double *matrix = NULL;
	double dsize;
	char out[50];
	strcpy(out,file);
	strcat(out,".temp");

//	double *zero=NULL;
//	int zero_size;
//	read_matrixB2(&zero, &zero_size, file);
//	show_matrix(zero, zero_size, file);

	if( !(f_file = fopen(file, "rb") ) )
	{
		printf("Msg(swapSPmatrix): I'm sorry, unable to open %s file\nForcing exit\n",file);
		exit(1);
	}
	fread(&dsize,sizeof(double),1,f_file); // getting size
	size = (long int) dsize; // casting double into integer
	if(debug)
	{
		printf("\nMsg(swapSPmatrix): File %s contains a %ld sized matrix\n",file,size);
		fflush(stdout);
	}

	// Some checking...
	if(size > mem)
	{
		printf("Msg(swapSPmatrix): Sorry, I can't speed up a %ld sized matrix and %ld elements memory.\nForcing exit\n",size,mem);
		exit(1);
	}

	if( !(f_out = fopen(out, "wb") ) )
	{
		printf("Msg(swapSPmatrix): I'm sorry, unable to open %s file\nForcing exit\n","out.bin");
		exit(1);
	}
	fwrite(&dsize,sizeof(double),1,f_out);  // write output matrix size

	// Stripe memory allocation
	if( !(matrix = (double *) malloc( sizeof( double ) * mem ) ) )
	{
		printf("Msg(swapSPmatrix): Sorry, unable to allocate %ld bytes\nForcing exit\n",sizeof(double) * stripe);
		exit(1);
	}

	// Loop stripes
	r = size; // remaining rows...
	n = 0; // stripe start
	m = size;

		// Loading stripe... (from disk)
		// (screening the whole Hessian matrix as it was computed and dumped into disk...)

		while( n < size )
		{
			fseek( f_file, 1*sizeof(double), SEEK_SET ); // positioning at file start
//			// Variable stripes...
//			if(r*(r+1)/2 > mem) // then stripes...
//				p = ( (2*r+1) - sqrt(4*r*r+1+4*r-8*block) ) / 2; // number of rows given "block" and "size"
//			else // last loop and exit...
//				p = r;

			p = mem / (size-n); // rectangular stripe ... (variable "p")
			if(p+n > size)
				p = size-n;

			if(debug)
				fprintf(stderr,"Loop %ld: p=%d rows stripe (for %ld elements block and r=%d rows matrix)\n",n,p,block,r);
			if(debug)
				fprintf(stderr,"Reading Stripe: rows= %ld-%ld   cols= %ld-%ld\n",n,n+p,0,size-n);

			// Read stripe...
			d = 0;
			c = 0; // stripe-matrix index counter
			for(y = 0; y < (size-n); y++) // (n+p "full" data chunks) + (p-1 smaller chunks)
			{
				fseek( f_file, n*sizeof(double), SEEK_CUR ); // positioning at stripe start
				if(y > (size-n-p))
					d++;
				for(x = 0; x < p-d; x++ )
				{
//					fread(&matrix[c],sizeof(double),1,f_file); // placing elements into memory
					fread(&matrix[y*p+x],sizeof(double),1,f_file); // placing elements into memory
					if(debug)
						fprintf(stderr,"%f ",matrix[y*p+x]);
//					c++;
				}
				if(debug)
					fprintf(stderr,"\n");
				next = size-(p+n+y);
				if(next > 0)
				{
					fseek( f_file, next*sizeof(double), SEEK_CUR ); // positioning at next column start
					if(debug)
						fprintf(stderr,"next= %ld\n",next);
				}
			}

			// Write stripe...
			if(debug)
				fprintf(stderr,"Writing matrix into file\n");
			for(x = 0; x < p; x++ )
			{
				for(y = size-n-x-1; y >= 0; y--) // <-- Reversing column order !!!
				{
					if(debug)
						fprintf(stderr,"%f ",matrix[y*p + x]);
					fwrite(&matrix[y*p + x],sizeof(double),1,f_out);
				}
				if(debug)
					fprintf(stderr,"\n");
			}
			n+=p;
		}

	// Closing files...
	fclose(f_file);
	if(debug)
		fprintf(stderr,"File %s closed\n",file);
	fclose(f_out);
	if(debug)
		fprintf(stderr,"File %s closed\n",out);

//	double *first=NULL,*second=NULL;
//	int first_size,second_size;
//	read_matrixB2(&first, &first_size, file);
//	show_matrix(first, first_size, file);
//	read_matrixB(&second, &second_size, out);
//	show_matrix(second, second_size, out);

	// Deleting input file
	if(debug)
		fprintf(stderr,"Deleting %s ... ",file);
	if(remove(file) == 0)
	{
		if(debug)
			fprintf(stderr,"File %s deleted successfully.\n",file);
	}
	else
		if(debug)
			fprintf(stderr, "Error deleting the file %s.\n",file);

	// Renaming output file
	if(debug)
		fprintf(stderr,"Renaming %s to %s ... ",out,file);
	if(rename(out,file) == 0)
	{
		if(debug)
			fprintf(stderr,"Success!\n");
	}
	else
		if(debug)
			fprintf(stderr, "Error renaming %s.\n", out);

}


// Computing the Radius of Gyration (Rg) of a Macromolecule.
// See: Seeliger and de Groot. PLOS (2010).
float radius_gyration(Macromolecule *mol)
{
	pdbIter *iter;
	Atom *atom;
	Tcoor pos;
	double mtot = 0.0;
	double mta = 0.0;
	double com[3];
	double radg[3];
	double radg2=0.0;
	iter = new pdbIter( mol );

	// Computing CoM
	// Shifting coordinates to their Center of Mass (Computing the CoM)
	// The structure should be centered on its CoM in order
	// to compute Kinetic-energy and Hessian Matrices.
	com[0] = com[1] = com[2] = 0.0;
	for ( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() )  // screens all-atoms
	{
		atom = iter->get_atom();
		atom->getPosition(pos);
		mta = atom->getPdbocc(); // Load mass from occupancies...
		mtot += mta;
		// Sum(mass*coord) before putting the CoM at 0
		com[0] += mta * pos[0];
		com[1] += mta * pos[1];
		com[2] += mta * pos[2];
	}
	com[0] /= mtot;
	com[1] /= mtot;
	com[2] /= mtot;
//	printf( "Msg(radius_gyration)> Mass %12.8f Center %.8f %.8f %.8f --> to 0,0,0\n", mtot, com[0], com[1], com[2] );

	// Computing Rg
	radg[0] = radg[1] = radg[2] = 0.0;
	for( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() )
	{
		atom = iter->get_atom();
		atom->getPosition(pos);
		mta = atom->getPdbocc(); // Load mass from occupancies...
		radg[0] = pos[0]-com[0]; // position respecting CoM
		radg[1] = pos[1]-com[1];
		radg[2] = pos[2]-com[2];
		radg[0] *= radg[0]; // square
		radg[1] *= radg[1];
		radg[2] *= radg[2];
		radg2 += (radg[0] + radg[1] + radg[2]) * mta;
	}
	delete iter;

	radg2 = sqrt( radg2/mtot );
	return((float) radg2);
}

// Length of one vector
double length(double *v1)
{
	return ( sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]) );
}

// Dot product between two vectors
double dotp(double *v1, double *v2)
{
	return ( v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] );
}
