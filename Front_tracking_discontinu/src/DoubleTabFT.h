/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
//////////////////////////////////////////////////////////////////////////////
//
// File:        DoubleTabFT.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////
#ifndef DoubleTabFT_included
#define DoubleTabFT_included
#include <DoubleTab.h>
#include <IntTab.h>
class DoubleTabFT : public DoubleTab
{
public:
  DoubleTabFT()
  {
    set_smart_resize(1);
  }
  DoubleTabFT(const DoubleTabFT&) = default;
  DoubleTabFT(int i, int j)
  {
    set_smart_resize(1);
    resize(i, j);
  }
  DoubleTabFT(const DoubleTab& x)
  {
    set_smart_resize(1);
    DoubleTab::operator=(x);
  }
  DoubleTabFT& operator=(const DoubleTab& x)
  {
    DoubleTab::operator=(x);
    return *this;
  }
  DoubleTabFT& operator=(const DoubleTabFT& x)
  {
    DoubleTab::operator=(x);
    return *this;
  }
  DoubleTabFT& operator=(const double x)
  {
    DoubleTab::operator=(x);
    return *this;
  }
};

class IntTabFT : public IntTab
{
public:
  IntTabFT()
  {
    set_smart_resize(1);
  }
  IntTabFT(const IntTabFT&) = default;
  IntTabFT(int i, int j)
  {
    set_smart_resize(1);
    resize(i, j);
  }
  IntTabFT& operator=(const IntTab& x)
  {
    IntTab::operator=(x);
    return *this;
  }
  IntTabFT& operator=(const IntTabFT& x)
  {
    IntTab::operator=(x);
    return *this;
  }
  IntTabFT& operator=(const int x)
  {
    IntTab::operator=(x);
    return *this;
  }
};

class ArrOfDoubleFT : public ArrOfDouble
{
public:
  ArrOfDoubleFT()
  {
    set_smart_resize(1);
  }
  ArrOfDoubleFT(int n)
  {
    set_smart_resize(1);
    resize_array(n);
  }
  ArrOfDoubleFT(const ArrOfDouble& x)
  {
    set_smart_resize(1);
    ArrOfDouble::operator=(x);
  }
  // on n 'appelle pas le constructeur par copie de ArrofDouble
  ArrOfDoubleFT(const ArrOfDoubleFT& x) : ArrOfDouble()
  {
    set_smart_resize(1);
    ArrOfDouble::operator=(x);
  }
  ArrOfDoubleFT& operator=(double i)
  {
    ArrOfDouble::operator=(i);
    return *this;
  }
};

class ArrOfIntFT : public ArrOfInt
{
public:
  ArrOfIntFT()
  {
    set_smart_resize(1);
  }
  ArrOfIntFT(const ArrOfIntFT&) = default;
  ArrOfIntFT(int n)
  {
    set_smart_resize(1);
    resize_array(n);
  }
  ArrOfIntFT& operator=(int i)
  {
    ArrOfInt::operator=(i);
    return *this;
  }
  ArrOfIntFT& operator=(const ArrOfInt& x)
  {
    ArrOfInt::operator=(x);
    return *this;
  }
  ArrOfIntFT& operator=(const ArrOfIntFT& x)
  {
    ArrOfInt::operator=(x);
    return *this;
  }
};
#endif
