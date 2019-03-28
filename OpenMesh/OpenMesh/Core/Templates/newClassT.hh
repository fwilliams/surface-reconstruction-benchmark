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
//  CLASS newClass
//
//=============================================================================
#ifndef DOXY_IGNORE_THIS
#ifndef OPENMESH_NEWCLASST_HH
#define OPENMESH_NEWCLASST_HH


//== INCLUDES =================================================================


//== FORWARDDECLARATIONS ======================================================


//== NAMESPACES ===============================================================

namespace OpenMesh {


//== CLASS DEFINITION =========================================================



	      
/** \class newClassT newClassT.hh <OpenMesh/.../newClassT.hh>

    Brief Description.
  
    A more elaborate description follows.
*/

template <>
class newClassT
{
public:
   
  /// Default constructor
  newClassT() {}
 
  /// Destructor
  ~newClassT() {}

  
private:

  /// Copy constructor (not used)
  newClassT(const newClassT& _rhs);

  /// Assignment operator (not used)
  newClassT& operator=(const newClassT& _rhs);
  
};


//=============================================================================
} // namespace OpenMesh
//=============================================================================
#if defined(OM_INCLUDE_TEMPLATES) && !defined(OPENMESH_NEWCLASS_C)
#define OPENMESH_NEWCLASS_TEMPLATES
#include "newClass.cc"
#endif
//=============================================================================
#endif // OPENMESH_NEWCLASST_HH defined
#endif // DOXY_IGNORE_THIS
//=============================================================================
