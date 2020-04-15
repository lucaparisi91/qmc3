

sinOrbital::sinOrbital(int n1,int n2,int n3,double lBox_) : lBox(lBox_) , ns{n1,n2,n3}
{
    k[0]=n1*2*M_PI/lBox;
    k[1]=n2*2*M_PI/lBox;
    k[2]=n3*2*M_PI/lBox;
    
    int isCos=1;
    
    if (abs(n1)>0 )
      {
	isCos=n1/std::abs(n1);
      }
    else
      {
	if (abs(n2)>0)
	  {
	    isCos=n2/std::abs(n2);
	  }
	else
	  {
	    if (abs(n3)>0)
	      {
		isCos=n3/std::abs(n3);
	      }
	  }
	  
      }
	
    if (isCos==1)
      {
	delta=M_PI/2.;
      }
    else
      {
	delta=0;
      }
   
  }


storeMatrixOrbitals(const state_t & states,matrix_t & matrix,std::vector<orbital_t> orbitals)
{
  const int D = states.dimensions()[1];
  
  for (const auto & orbital : orbitals)
    {
      orbital(states)
    }
}
