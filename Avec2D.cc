/*********************************************************************
Copyright 2000 - 2006 Bruno Wittmer
          2006 Michael Henke
This file is part of the Avec package.

Avec is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Avec is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Avec; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**********************************************************************/

#include "Avec2D.h"
#include <iostream>
#include <iomanip>
#include <iterator>


ClassImp(Avec2D)

// Transform Avec into Avec2D with given number of columns
Avec2D::Avec2D(const Avec & source, size_type cols)
{
  size_type rows = source.size()/cols;
  //std::cout << "Filling Avec2D with " 
  //	    << cols << " cols and "
  //	    << rows << " rows" << std::endl;
  for(size_type f = 0; f < cols; f++){
    push_back(Avec(rows));
    for(size_type g = 0; g < rows ; g++){
      size_type index=g+f*rows;
      (*this)[f][g]=source[index];
    }
  }
}

// Avec2D& Avec2D::operator=(const Avec2D& in)
// {
//   *this = in;
//   return *this;
// }


//! Print contents to a stream
void Avec2D::print(std::ostream & out, const char* sep) const
{
  out << std::setprecision(8);
  std::copy(begin(),end(),std::ostream_iterator<Avec>(out,sep));
}

std::ostream& operator<< (std::ostream& o, const Avec2D& v)
{
  (vtrans(v)).print(o);
  return o;
}

//! Concatenate vectors
Avec2D operator&(const Avec2D& vec1, const Avec2D& vec2)
{
  Avec2D result;
  if(vec1.size() != vec2.size()){
    std::cerr << "Cannot concatenate 2D Avecs wit different number of cols\n";
    return result;
  }
  Avec2D::const_iterator i1=vec1.begin();
  Avec2D::const_iterator i2=vec2.begin();
  for(;i1 != vec1.end();i1++,i2++)
    result.push_back(*i1 & *i2);
  return result;
}

//! Write vector to a file
void vecwrite(const Avec2D &vec, const std::string &filename)
{
  if(Avec::VERBOSE_FLAG) std::cout << "Writing 2D vector to file " << filename << std::endl;
  std::fstream file(filename.c_str(), std::ios::out);
  if(!file)return;
  file << vec;
  file.close();
}

//! Read vector from a file
int vecread(Avec2D &vec, const std::string &filename)
{
  std::fstream infile(filename.c_str(), std::ios::in);
  if(!infile){
    std::cerr << "Could not open file " << filename << std::endl;
    return -1;
  }
  
  if(vec.size() != 0){
    // Size is predefined, read only the amount of elements needed
    std::vector<Avec*> list;
    for(std::vector<Avec>::iterator i= vec.begin(); i != vec.end(); i++)
      list.push_back(&(*i));
    //{
    //  for(Avec::iterator f = i->begin(); f != i->end(); f++){
    //	if(! (infile>> *f))break;
    //  }
    //}
    vecread(list,filename);
  }
  else{
    // Vector is empty, assume that file is formatted in a table
    //vec.clear();

    dattype b=0;
    int newline=0;

    while( infile >> b){
      vec.push_back(Avec(1,b));
      char a;
      while (infile.get(a)){
	if( !isspace(a)) {
	  //std::cout << "\t";
	  infile.putback(a);
	  break;
	}
	if(a == '\n'){
	  newline=1;
	}
      }
      if(newline)break;
    }
    //std::cout << endl;
    while(infile){
      for(Avec2D::iterator i=vec.begin();i != vec.end(); i++){
	if(infile >> std::skipws >> b)
	  i->push_back(b);
	else
	  break;
      }
    }
  }
    
  if(Avec::VERBOSE_FLAG){
    std::cout << "Read in 2D vector with " << vec.size() << " columns";
    if(vec.size() != 0)
      std::cout << " and " << vec[0].size() << " rows";
    std::cout << std::endl;
  }
  infile.close();
  return 0;
}

//! Maxima of all columns
Avec vmaxc(const Avec2D& v1)
{
  Avec::size_type inp_col = v1.size();

  Avec v2(inp_col,0.0);
  for( Avec::size_type i = 0; i < inp_col; i++){
      v2[i] = vmax(v1[i]);
  }
  return v2;
}

//! Maxima of all rows
Avec vmaxr(const Avec2D& v1)
{
  return vmaxc(vtrans(v1));
}

//! Maximum of all elements
dattype vmax(const Avec2D &vec)
{
  return vmax(vmaxc(vec));
}

//! Minima of all columns
Avec vminc(const Avec2D& v1)
{
  Avec::size_type inp_col = v1.size();

  Avec v2(inp_col,0.0);
  for( Avec::size_type i = 0; i < inp_col; i++){
      v2[i] = vmin(v1[i]);
  }
  return v2;
}

//! Minima of all rows
Avec vminr(const Avec2D& v1)
{
  return vminc(vtrans(v1));
}

//! Minimum of all elements
dattype vmin(const Avec2D &vec)
{
  return vmin(vminc(vec));
}

//! Sum of all elements
dattype vsum(const Avec2D& v1)
{
  return(vsum(vsumc(v1)));
}

//! function to sum up an Avec2D in one dimension. 
/*! In this case the single Avecs(="columns") are summed up 
    one by one. As output one gets an Avec with the length 
    equals to the number of Avecs, in each position the total 
    sum of the Avec in this position in the input Avec2D
*/
Avec vsumc(const Avec2D& v1)
{    
  Avec::size_type inp_col = v1.size();

  Avec v2(inp_col,0.0);
  for( Avec::size_type i = 0; i < inp_col; i++)
    {
      v2[i] = vsum(v1[i]);
    }
  return v2;
}

//! function to sum up an Avec2D in one dimension. 
/*  In this case the "rows" are summed up one by one. 
    As output on gets an Avec with the length equals to the 
    length of one Avec, in each position the total sum of the 
    doubles in this position in the input Avec2D
*/
Avec vsumr(const Avec2D& v1)
{
  return vsumc(vtrans(v1));
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//##############################################################################################
//# function to "transpone" an Avec2D                                                          #
//##############################################################################################

Avec2D vtrans(const Avec2D& v1)
{
  Avec::size_type inp_row = v1.size();
  Avec::size_type inp_col = v1[0].size();

  Avec2D v2(inp_col,inp_row);
  for( Avec::size_type col = 0; col < inp_col; col++){
    for( Avec::size_type row = 0; row < inp_row; row++){
      v2[col][row] = v1[row][col];
    }
  }
  return v2;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Resize the Avec2D in both dimensions
Avec2D avec_resize(const Avec2D& v, Avec2D::size_type nsizex, Avec::size_type nsizey)
{
  Avec2D result;
  if(v.empty() || nsizex == 0)return result;
  if(nsizey == 0)nsizey=v[0].size();

  if(Avec::VERBOSE_FLAG) std::cout << "New size: " << nsizex << " x " << nsizey << std::endl;

  result = Avec2D(v.size(), nsizey);

  for(Avec2D::size_type i=0; i!=v.size();i++){
    result[i]=avec_resize(v[i],nsizey);
  }

  result=vtrans(result);
  Avec2D result2(nsizey, nsizex);

  for(Avec2D::size_type i=0; i!=result.size();i++){
    result2[i]=avec_resize(result[i],nsizex);
  }

  return vtrans(result2);
}

// Resize the Avec2D in both dimensions by factor
//Avec2D rescale_size(const Avec2D& v, double factorx, double factory)
//{
//  Avec2D result;
//  if(v.empty())return result;

//  nx=(int)((double)v.size()*factor);

//  return result;
//}

/////////////////////////////
// Arithmetic with scalars //
/////////////////////////////

Avec2D Avec2D::operator+(double offset) const
{
  Avec2D result(size());

  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    *f=*i+offset;
  return result;
}

Avec2D Avec2D::operator-(double offset) const
{
  Avec2D result(size());

  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    *f=*i-offset;
  return result;
}

Avec2D Avec2D::operator*(dattype factor) const
{
  Avec2D result(size());

  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    *f=*i*factor;
  return result;
}

Avec2D Avec2D::operator/(double factor) const
{
  Avec2D result(size());
  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    *f=*i/factor;
  return result;
}

void Avec2D::operator+=(double offset)
{
  for(iterator i=begin();i!=end();i++)
    *i+=offset;
}

void Avec2D::operator-=(double offset)
{
  for(iterator i=begin();i!=end();i++)
    *i-=offset;
}

void Avec2D::operator*=(double factor)
{
  for(iterator i=begin();i!=end();i++)
    *i*=factor;
}

void Avec2D::operator/=(double factor)
{
  for(iterator i=begin();i!=end();i++)
    *i/=factor;
}

Avec2D Avec2D::operator-() const
{
  return *this*-1;
}

////////////////////////////////
// Arithmetic with 2D-vectors //
////////////////////////////////

Avec2D Avec2D::operator+(const Avec2D& source) const
{
  size_t length=std::min(size(),source.size());
  Avec2D result(length);

  const_iterator f=begin();
  const_iterator g=source.begin();
  for(iterator i=result.begin();i!=result.end();i++,f++,g++)
    *i=*f+*g;
  return result;
}

Avec2D Avec2D::operator-(const Avec2D& source) const
{
  size_t length=std::min(size(),source.size());
  Avec2D result(length);

  const_iterator f=begin();
  const_iterator g=source.begin();
  for(iterator i=result.begin();i!=result.end();i++,f++,g++)
    *i=*f-*g;
  return result;
}

Avec2D Avec2D::operator*(const Avec2D& source) const
{
  size_t length=std::min(size(),source.size());
  Avec2D result(length);

  const_iterator f=begin();
  const_iterator g=source.begin();
  for(iterator i=result.begin();i!=result.end();i++,f++,g++)
    *i=(*f)*(*g);
  return result;
}

Avec2D Avec2D::operator/(const Avec2D& source) const
{
  size_t length=std::min(size(),source.size());
  Avec2D result(length);

  const_iterator f=begin();
  const_iterator g=source.begin();
  for(iterator i=result.begin();i!=result.end();i++,f++,g++)
    *i=(*f)/(*g);
  return result;
}

void Avec2D::operator+=(const Avec2D& source)
{
  size_t length=std::min(size(),source.size());

  const_iterator g=source.begin();
  for(iterator i=begin();i!=begin()+length;i++,g++)
    *i+=*g;
}

void Avec2D::operator-=(const Avec2D& source)
{
  size_t length=std::min(size(),source.size());

  const_iterator g=source.begin();
  for(iterator i=begin();i!=begin()+length;i++,g++)
    *i-=*g;
}

void Avec2D::operator*=(const Avec2D& source)
{
  size_t length=std::min(size(),source.size());

  const_iterator g=source.begin();
  for(iterator i=begin();i!=begin()+length;i++,g++)
    *i*=*g;
}

void Avec2D::operator/=(const Avec2D& source)
{
  size_t length=std::min(size(),source.size());

  const_iterator g=source.begin();
  for(iterator i=begin();i!=begin()+length;i++,g++)
    *i/=*g;
}

////////////////////////////////
// Arithmetic with 1D-vectors //
////////////////////////////////

Avec2D Avec2D::operator+(const Avec& vec) const
{
  Avec2D result(size());

  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    *f=*i+vec;
  return result;
}

Avec2D Avec2D::operator-(const Avec& vec) const
{
  Avec2D result(size());

  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    *f=*i-vec;
  return result;
}

Avec2D Avec2D::operator*(const Avec& vec) const
{
  Avec2D result(size());

  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    *f=*i*vec;
  return result;
}

Avec2D Avec2D::operator/(const Avec& vec) const
{
  Avec2D result(size());
  iterator f=result.begin();
  for(const_iterator i=begin();i!=end();i++,f++)
    *f=*i/vec;
  return result;
}

void Avec2D::operator+=(const Avec& vec)
{
  for(iterator i=begin();i!=end();i++)
    *i+=vec;
}

void Avec2D::operator-=(const Avec& vec)
{
  for(iterator i=begin();i!=end();i++)
    *i-=vec;
}

void Avec2D::operator*=(const Avec& vec)
{
  for(iterator i=begin();i!=end();i++)
    *i*=vec;
}

void Avec2D::operator/=(const Avec& vec)
{
  for(iterator i=begin();i!=end();i++)
    *i/=vec;
}


// Comparisons
Avec2D Avec2D::operator<(const dattype value) const
{
  size_type cols=size();
  Avec2D result(cols);

  for(size_type i=0 ;i!=cols;i++)
    result[i]=(*this)[i]<value;
  return result;
}

Avec2D Avec2D::operator>(const dattype value) const
{
  size_type cols=size();
  Avec2D result(cols);

  for(size_type i=0 ;i!=cols;i++)
    result[i]=(*this)[i]>value;
  return result;
}

Avec2D Avec2D::operator!=(const dattype value) const
{
  size_type cols=size();
  Avec2D result(cols);

  for(size_type i=0 ;i!=cols;i++)
    result[i]=(*this)[i]!=value;
  return result;
}

Avec2D Avec2D::operator&&(const Avec2D& source) const
{
  size_type cols=size();
  Avec2D result(cols);

  for(size_type i=0 ;i!=cols;i++)
    result[i]=(*this)[i]&&source[i];
  return result;
}

Avec2D Avec2D::operator||(const Avec2D& source) const
{
  size_type cols=size();
  Avec2D result(cols);

  for(size_type i=0 ;i!=cols;i++)
    result[i]=(*this)[i]||source[i];
  return result;
}


// math.h functions
Avec2D acos(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(acos(*i));
  return result;
}

Avec2D asin(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(asin(*i));
  return result;
}

Avec2D atan(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(atan(*i));
  return result;
}

Avec2D atan2(const Avec2D &y, const Avec2D &x)
{
  Avec2D result;
  Avec2D::const_iterator i,f;
  for(i=y.begin(),f=x.begin();i!=y.end() && f!=x.end();i++,f++)
    result.push_back(atan(*i));
  return result;
}

Avec2D ceil(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(ceil(*i));
  return result;
}

Avec2D cos(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(cos(*i));
  return result;
}

Avec2D cosh(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(cosh(*i));
  return result;
}

Avec2D exp(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(exp(*i));
  return result;
}

WINAVEC_API Avec2D fabs(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(fabs(*i));
  return result;
}

Avec2D floor(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(floor(*i));
  return result;
}

Avec2D fmod(const Avec2D &vec, dattype y)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(fmod(*i, y));
  return result;
}

Avec2D log(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(log(*i));
  return result;
}


Avec2D log10(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(log10(*i));
  return result;
}


Avec2D pow(const Avec2D &vec, dattype y)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(pow(*i, y));
  return result;
}

Avec2D sin(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(sin(*i));
  return result;
}

Avec2D sinh(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(sinh(*i));
  return result;
}


Avec2D sqrt(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(sqrt(*i));
  return result;
}

Avec2D tan(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(tan(*i));
  return result;
}

Avec2D tanh(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(tanh(*i));
  return result;
}


//////////////////////////
//    Random Numbers    //
//////////////////////////

//! Random numbers from 0 to RAND_MAX
Avec2D rand(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(rand(*i));
  return result;
}

Avec2D urand(const Avec2D &vec)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(urand(*i));
  return result;
}

Avec2D grand(const Avec2D &vec, double sigma)
{
  Avec2D result;
  for(Avec2D::const_iterator i=vec.begin();i!=vec.end();i++)
    result.push_back(grand(*i, sigma));
  return result;
}


void avec_remove_nan(Avec2D& data, dattype replacement)
{
  for(Avec2D::iterator i = data.begin(); i != data.end(); i++)
    avec_remove_nan(*i, replacement);
}


/////////////////////////////////////////////////////////////////////////

#ifdef __AVECROOT__
// Root-specific stuff
#include "TH2.h"
#include "TDirectory.h"
#include "TVirtualPad.h"
#include "TMultiGraph.h"
#include "TFile.h"

// Overloaded TObject methods
void Avec2D::Draw(Option_t* option)
{
  avec_draw(*this);
}

void Avec2D::Browse(TBrowser* b)
{
  gPad->Clear();
  avec_draw(*this);
  gPad->Update();
}

Avec2D empty_avec2d;

//! Retrieve an Avec2D from a root file through TFile reference
Avec2D& avec2d_get(const std::string& name, TFile& file)
{
  Avec2D* Avec_ptr = 0;
  file.GetObject(name.c_str(),Avec_ptr);
  if(!Avec_ptr){
    std::cerr << "Did not find " << name << std::endl;
    return empty_avec2d;
  }
  return *Avec_ptr;
}

//! Retrieve an Avec2D from a root file through file name
Avec2D& avec2d_get(const std::string& vecname, const std::string& filename)
{
  TFile f(filename.c_str(),"READ");
  if(!f.IsOpen()){
    std::cerr << "Could not open file " << filename << std::endl;
    return empty_avec2d;
  }
  return avec2d_get(vecname, f);
}

// Draw columns of Avec2D against one Avec (must be of same size)
WINAVEC_API TMultiGraph *avec_draw(const Avec& x, const Avec2D& y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  if((x.size() == 0) || (y.size() == 0)) return 0;

  // Put all columns into Avec list
  std::vector<const Avec*> list;
  for(Avec2D::size_type i=0; i!=y.size(); i++){
    // Check if sets are of same size      
    if (x.size() != y[i].size()){
      std::cerr << "avec_draw: Column " << i << " has not same size as x-vector!" << std::endl;
    }
    else list.push_back(&(y[i]));
  }

  // reuse function for drawing Avec list
  return avec_draw(x, list, title, xTitle, yTitle, options, draw, name);
}

// Draw columns of Avec2D with errors against one Avec (must be of same size)
WINAVEC_API TMultiGraph *avec_draw(const Avec& x, const Avec2D& y, const Avec2D& ey, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  if((x.size() == 0) || (y.size() == 0) || (ey.size() == 0)) return 0;

  // Put all columns into Avec list
  std::vector<const Avec*> list;
  std::vector<const Avec*> err_list;
  for(Avec2D::size_type i=0; i!=y.size(); i++){
    // Check if sets are of same size      
    if (x.size() != y[i].size() || x.size() != ey[i].size()){
      std::cerr << "avec_draw: Column " << i << " has not same size as x-vector!" << std::endl;
    }
    else{
    	 list.push_back(&(y[i]));
    	 err_list.push_back(&(ey[i]));
    }
  }

  // reuse function for drawing Avec list
  return avec_draw(x, list, err_list, title, xTitle, yTitle, options, draw, name);
}


// Draw columns of two Avec2D against each other
WINAVEC_API TMultiGraph *avec_draw(const Avec2D& x, const Avec2D& y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  if((x.size() == 0) || (y.size() == 0)) return 0;

  if(x.size() != y.size()){
    std::cerr << "Error in avec_draw : The two Avec2D are not of the same size" << std::endl;
    return 0;
  }

  // Put all columns into Avec list
  std::vector<const Avec*> list_x;
  std::vector<const Avec*> list_y;
  for(Avec2D::size_type i=0; i!=y.size(); i++){
    list_x.push_back(&(x[i]));
    list_y.push_back(&(y[i]));
  }

  // reuse function for drawing Avec list
  return avec_draw(list_x, list_y, title, xTitle, yTitle, options, draw, name);
}

// Draw columns of two Avec2D against each other
WINAVEC_API TMultiGraph *avec_draw(const Avec2D& x, const Avec2D& y, const Avec2D& err_y, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, const std::string & name)
{
  if((x.size() == 0) || (y.size() == 0) || (err_y.size() == 0)) return 0;

  if(x.size() != y.size() || y.size() != err_y.size()){
    std::cerr << "Error in avec_draw : The two Avec2D or the Error Avec2D are not of the same size" << std::endl;
    return 0;
  }

  // Put all columns into Avec list
  std::vector<Avec> list_x;
  std::vector<Avec> list_y;
  //std::vector<Avec> list_ex;
  std::vector<Avec> list_ey;
  for(Avec2D::size_type i=0; i!=y.size(); i++){
    list_x.push_back((x[i]));
    //list_ex.push_back( new Avec(x[i].size()));
    list_y.push_back((y[i]));
    list_ey.push_back((err_y[i]));
  }

  // reuse function for drawing Avec list
  return avec_draw(list_x, list_y, list_ey, title, xTitle, yTitle, options, draw, name);
}



WINAVEC_API TH2D *avec_draw(const Avec2D& z, const std::string& title, const std::string & xTitle, const std::string & yTitle, const std::string& options, int draw, std::string name)
{
  if(z.empty())return 0;
  if(name == "avec_draw" && title != "") name=title;

  if(Avec::VERBOSE_FLAG) std::cout << "name.c_str() is:    " << name.c_str() << std::endl;
      
  int nxbins = z.size();
  int nybins = z[0].size();


  std::string histo_name = GetValidObjectName<TH2D>(name);
  if(Avec::VERBOSE_FLAG) std::cout << "TH2D name is:    " << histo_name.c_str() << std::endl;

  TH2D* h1 = new TH2D( histo_name.c_str() , title.c_str() , nxbins , 0.5 , nxbins + 0.5 , nybins , 0.5 , nybins + 0.5 );
 
  for(int i=0; i<nxbins; i++){
    for(int f=0; f<nybins; f++){   
      h1->SetBinContent( i+1 , f+1, z[i][f]);
    }
  }

  h1->GetXaxis()->SetTitle(xTitle.c_str());
  h1->GetYaxis()->SetTitle(yTitle.c_str());
  //std::cout << "DrawOption: " << options << std::endl;
  h1->SetDrawOption(options.c_str());
  h1->Draw(options.c_str());
 
  return h1; 

}

//! Fill all values of one Avec2D into a histogram weighting with a second Avec2D
/*! If w is an empty Avec2D default weight is 1
*/
void hfill(const Avec2D& v, TH1& h, const Avec2D& w)
{
  int w_flag = (w.size()!=0);
  if(v.size() != w.size() &&   w_flag){
    std::cerr << "Error in hfill(const Avec2D&, TH1&, const Avec2D&): The two Avec2D should have the same size!\n";
    return;
  }
  for(Avec2D::size_type i=0; i<v.size(); i++){
    Avec weight;
    if(w_flag)weight=w[i];
    else weight=Avec(v[i].size(),1);
    hfill(v[i], h, weight);
  }
}

// Fill the conetents of an Avec into a histogram and plot it
TH1* avec_plot(const Avec2D& v, int nbins, const std::string& title, const std::string& name)
{
  // Create the limits for the x-axis
  double xmax=vmax(v);
  double xmin=vmin(v);
  double diff=xmax-xmin;
  xmax+= 0.1*diff;
  xmin-= 0.1*diff;

  // Check if TH1 object with this name already exists, if yes remove it
  TH1* h_test = (TH1*)(gDirectory->FindObject( name.c_str()));
  if( h_test != 0 ){
    h_test->Delete();
    if(Avec::VERBOSE_FLAG) std::cout << "existing TH1 deleted!" << std::endl;
  }

  TH1 * h=new TH1D(name.c_str(),title.c_str(),nbins,xmin,xmax);

  hfill(v,*h);

  h->Draw();
  return h;
}



#endif
