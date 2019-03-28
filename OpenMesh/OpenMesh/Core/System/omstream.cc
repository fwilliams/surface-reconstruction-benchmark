/*===========================================================================*\
 *                                                                           *
 *                               OpenMesh                                    *
 *      Copyright (C) 2001-2003 by Computer Graphics Group, RWTH Aachen      *
 *                           www.openmesh.org                                *
 *                                                                           *
 *---------------------------------------------------------------------------* 
 *                                                                           *
 *                                License                                    *
 *                                                                           *
 *  This library is free software; you can redistribute it and/or modify it  *
 *  under the terms of the GNU Lesser General Public License as published    *
 *  by the Free Software Foundation, version 2.1.                            *
 *                                                                           *
 *  This library is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        *
 *  Lesser General Public License for more details.                          *
 *                                                                           *
 *  You should have received a copy of the GNU Lesser General Public         *
 *  License along with this library; if not, write to the Free Software      *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                *
 *                                                                           *
\*===========================================================================*/

//=============================================================================
//
//  CLASS mostream - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================

#include <OpenMesh/Core/System/omstream.hh>


//== IMPLEMENTATION ========================================================== 


OpenMesh::mostream& omlog() 
{
  static bool initialized = false;
  static OpenMesh::mostream mystream;
  if (!initialized)
  {
    mystream.connect(std::clog);
#ifdef NDEBUG
    mystream.disable();
#endif
  }
  return mystream;
}


OpenMesh::mostream& omout() 
{
  static bool initialized = false;
  static OpenMesh::mostream mystream;
  if (!initialized) mystream.connect(std::cout);
  return mystream;
}


OpenMesh::mostream& omerr() 
{
  static bool initialized = false;
  static OpenMesh::mostream mystream;
  if (!initialized) mystream.connect(std::cerr);
  return mystream;
}


//=============================================================================
