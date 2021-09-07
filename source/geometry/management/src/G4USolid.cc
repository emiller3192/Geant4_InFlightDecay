//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id:$
// GEANT4 tag $Name:$
//
//
// G4USolid implementation
//
// --------------------------------------------------------------------

#include "G4USolid.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"
#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4VisExtent.hh"
#include "G4PhysicalConstants.hh"
#include "G4GeometryTolerance.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

G4USolid::G4USolid(const G4String& name, VUSolid* s) :
  G4VSolid(name), fShape(s), fRebuildPolyhedron(false), fPolyhedron(0)
{
}

G4USolid::G4USolid(__void__& a)
  : G4VSolid(a), fShape(0), fRebuildPolyhedron(false), fPolyhedron(0)
{
}

G4USolid::~G4USolid()
{
  delete fPolyhedron; fPolyhedron = 0;
}

G4bool G4USolid::operator==(const G4USolid& s) const
{
  return (this == &s) ? true : false;
}

EInside G4USolid::Inside(const G4ThreeVector& p) const
{
  UVector3 pt;
  VUSolid::EnumInside in_temp;
  EInside in = kOutside;
  pt.x() = p.x();
  pt.y() = p.y();
  pt.z() = p.z(); // better assign at construction

  in_temp = fShape->Inside(pt);

  if (in_temp == VUSolid::EnumInside::eSurface) return kSurface;
  if (in_temp == VUSolid::EnumInside::eInside) return kInside;

  return in;
}

G4ThreeVector G4USolid::SurfaceNormal(const G4ThreeVector& pt) const
{
  UVector3 p;
  p.x() = pt.x();
  p.y() = pt.y();
  p.z() = pt.z();
  UVector3 n;
  fShape->Normal(p, n);
  return G4ThreeVector(n.x(), n.y(), n.z());
}

G4double G4USolid::DistanceToIn(const G4ThreeVector& pt,
                                const G4ThreeVector& d) const
{
  UVector3 p;
  p.x() = pt.x();
  p.y() = pt.y();
  p.z() = pt.z(); // better assign at construction
  UVector3 v;
  v.x() = d.x();
  v.y() = d.y();
  v.z() = d.z(); // better assign at construction
  G4double dist = fShape->DistanceToIn(p, v);
  if (dist > kInfinity) return kInfinity;
//  return  (dist > halfTolerance) ? dist : 0.0;
  return dist;
}

G4double G4USolid::DistanceToIn(const G4ThreeVector& pt) const
{
  UVector3 p;
  p.x() = pt.x();
  p.y() = pt.y();
  p.z() = pt.z(); // better assign at construction
  G4double dist = fShape->SafetyFromOutside(p); // true?
  if (dist > kInfinity) return kInfinity;
//  return (dist > halfTolerance) ? dist : 0.0;
  return dist;
}

G4double G4USolid::DistanceToOut(const G4ThreeVector& pt,
                                 const G4ThreeVector& d,
                                 const G4bool calcNorm,
                                 G4bool* validNorm,
                                 G4ThreeVector* norm) const
{
  UVector3 p;
  p.x() = pt.x();
  p.y() = pt.y();
  p.z() = pt.z(); // better assign at construction
  UVector3 v;
  v.x() = d.x();
  v.y() = d.y();
  v.z() = d.z(); // better assign at construction
  UVector3 n;
  G4bool valid;
  G4double dist = fShape->DistanceToOut(p, v, n, valid); // should use local variable
  if(calcNorm)
  {
    if(valid){ *validNorm = true; }
    else { *validNorm = false; }
    if(*validNorm)  // *norm = n, but only after calcNorm check
    {
      norm->setX(n.x());
      norm->setY(n.y());
      norm->setZ(n.z());
    }
  }
  if (dist > kInfinity) return kInfinity;
//  return (dist > halfTolerance) ? dist : 0.0;
  return dist;
}

G4double G4USolid::DistanceToOut(const G4ThreeVector& pt) const
{
  UVector3 p;
  p.x() = pt.x();
  p.y() = pt.y();
  p.z() = pt.z(); // better assign at construction
  G4double dist = fShape->SafetyFromInside(p); // true?
//  return (dist > halfTolerance) ? dist : 0.0;
  return dist;
}

G4double G4USolid::GetCubicVolume()
{
  return fShape->Capacity();
}

G4double G4USolid::GetSurfaceArea()
{
  return fShape->SurfaceArea();
}

G4ThreeVector G4USolid::GetPointOnSurface() const
{
  UVector3 p;
  p = fShape->GetPointOnSurface();
  return G4ThreeVector(p.x(), p.y(), p.z());
}

G4bool G4USolid::CalculateExtent(const EAxis pAxis,
                                 const G4VoxelLimits& pVoxelLimit,
                                 const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const
{
  if (!pTransform.IsRotated())
  {
    VUSolid::EAxisType eAxis = VUSolid::eXaxis;
    G4double offset = pTransform.NetTranslation().x();
    if (pAxis == kYAxis)
    {
      eAxis = VUSolid::eYaxis;
      offset = pTransform.NetTranslation().y();
    }
    if (pAxis == kZAxis)
    {
      eAxis = VUSolid::eZaxis;
      offset = pTransform.NetTranslation().z();
    }
    fShape->ExtentAxis(eAxis, pMin, pMax);

    pMin += offset;
    pMax += offset;

    if (pVoxelLimit.IsLimited())
    {
      switch (pAxis)
      {
        case kXAxis:
          if ((pMin > pVoxelLimit.GetMaxXExtent() + kCarTolerance) ||
              (pMax < pVoxelLimit.GetMinXExtent() - kCarTolerance))
          {
            return false;
          }
          else
          {
            pMin = std::max(pMin, pVoxelLimit.GetMinXExtent());
            pMax = std::min(pMax, pVoxelLimit.GetMaxXExtent());
          }
          break;
        case kYAxis:
          if ((pMin > pVoxelLimit.GetMaxYExtent() + kCarTolerance) ||
              (pMax < pVoxelLimit.GetMinYExtent() - kCarTolerance))
          {
            return false;
          }
          else
          {
            pMin = std::max(pMin, pVoxelLimit.GetMinYExtent());
            pMax = std::min(pMax, pVoxelLimit.GetMaxYExtent());
          }
          break;
        case kZAxis:
          if ((pMin > pVoxelLimit.GetMaxZExtent() + kCarTolerance) ||
              (pMax < pVoxelLimit.GetMinZExtent() - kCarTolerance))
          {
            return false;
          }
          else
          {
            pMin = std::max(pMin, pVoxelLimit.GetMinZExtent());
            pMax = std::min(pMax, pVoxelLimit.GetMaxZExtent());
          }
          break;
        default:
          break;
      }
      pMin -= kCarTolerance ;
      pMax += kCarTolerance ;
    }
    return true;
  }
  else  // General rotated case - create and clip mesh to boundaries
  {
    // Rotate BoundingBox and Calculate Extent as for BREPS

    G4bool existsAfterClip = false ;
    G4ThreeVectorList* vertices ;

    pMin = +kInfinity ;
    pMax = -kInfinity ;

    // Calculate rotated vertex coordinates

    vertices = CreateRotatedVertices(pTransform) ;
    ClipCrossSection(vertices, 0, pVoxelLimit, pAxis, pMin, pMax) ;
    ClipCrossSection(vertices, 4, pVoxelLimit, pAxis, pMin, pMax) ;
    ClipBetweenSections(vertices, 0, pVoxelLimit, pAxis, pMin, pMax) ;

    if (pVoxelLimit.IsLimited(pAxis) == false)
    {
      if ((pMin != kInfinity) || (pMax != -kInfinity))
      {
        existsAfterClip = true ;

        // Add 2*tolerance to avoid precision troubles

        pMin -= kCarTolerance;
        pMax += kCarTolerance;
      }
    }
    else
    {
      G4ThreeVector clipCentre(
        (pVoxelLimit.GetMinXExtent() + pVoxelLimit.GetMaxXExtent()) * 0.5,
        (pVoxelLimit.GetMinYExtent() + pVoxelLimit.GetMaxYExtent()) * 0.5,
        (pVoxelLimit.GetMinZExtent() + pVoxelLimit.GetMaxZExtent()) * 0.5);

      if ((pMin != kInfinity) || (pMax != -kInfinity))
      {
        existsAfterClip = true ;


        // Check to see if endpoints are in the solid

        clipCentre(pAxis) = pVoxelLimit.GetMinExtent(pAxis);

        if (Inside(pTransform.Inverse().TransformPoint(clipCentre)) != kOutside)
        {
          pMin = pVoxelLimit.GetMinExtent(pAxis);
        }
        else
        {
          pMin -= kCarTolerance;
        }
        clipCentre(pAxis) = pVoxelLimit.GetMaxExtent(pAxis);

        if (Inside(pTransform.Inverse().TransformPoint(clipCentre)) != kOutside)
        {
          pMax = pVoxelLimit.GetMaxExtent(pAxis);
        }
        else
        {
          pMax += kCarTolerance;
        }
      }

      // Check for case where completely enveloping clipping volume
      // If point inside then we are confident that the solid completely
      // envelopes the clipping volume. Hence set min/max extents according
      // to clipping volume extents along the specified axis.

      else if (Inside(pTransform.Inverse().TransformPoint(clipCentre))
               != kOutside)
      {
        existsAfterClip = true ;
        pMin            = pVoxelLimit.GetMinExtent(pAxis) ;
        pMax            = pVoxelLimit.GetMaxExtent(pAxis) ;
      }
    }
    delete vertices;
    return existsAfterClip;
  }
}

void G4USolid::ComputeDimensions(G4VPVParameterisation*,
                                 const G4int,
                                 const G4VPhysicalVolume*)
{
    std::ostringstream message;
    message << "Illegal call to G4USolid::ComputeDimensions()" << G4endl
            << "Method not overloaded by derived class !";
    G4Exception("G4USolid::ComputeDimensions()", "GeomSolids0003",
                FatalException, message);
}

void G4USolid::DescribeYourselfTo(G4VGraphicsScene& scene) const
{
  scene.AddSolid(*this);
}

G4GeometryType G4USolid::GetEntityType() const
{

  G4String string = fShape->GetEntityType();
  return "G4" + string;
}

std::ostream& G4USolid::StreamInfo(std::ostream& os) const
{
  return fShape->StreamInfo(os);
}

G4USolid::G4USolid(const G4USolid& rhs)
  : G4VSolid(rhs), fRebuildPolyhedron(false), fPolyhedron(0)
{
  fShape = rhs.fShape->Clone();
}

G4USolid& G4USolid::operator=(const G4USolid& rhs)
{
  // Check assignment to self
  //
  if (this == &rhs)
  {
    return *this;
  }

  // Copy base class data
  //
  G4VSolid::operator=(rhs);

  // Copy data
  //
  fShape = rhs.fShape->Clone();
  fRebuildPolyhedron = false;
  delete fPolyhedron; fPolyhedron = 0;

  return *this;
}

G4VSolid* G4USolid::Clone() const
{
  std::ostringstream message;
  message << "Clone() method not implemented for type: "
          << GetEntityType() << "!" << G4endl
          << "Returning NULL pointer!";
  G4Exception("G4USolid::Clone()", "GeomSolids1001", JustWarning, message);
  return 0;
}

G4ThreeVectorList*
G4USolid::CreateRotatedVertices(const G4AffineTransform& pTransform) const
{
  G4double xMin, xMax, yMin, yMax, zMin, zMax;

  fShape->ExtentAxis(VUSolid::eXaxis, xMin, xMax);
  fShape->ExtentAxis(VUSolid::eYaxis, yMin, yMax);
  fShape->ExtentAxis(VUSolid::eZaxis, zMin, zMax);

  G4ThreeVectorList* vertices;
  vertices = new G4ThreeVectorList();

  if (vertices)
  {
    vertices->reserve(8);
    G4ThreeVector vertex0(xMin, yMin, zMin);
    G4ThreeVector vertex1(xMax, yMin, zMin);
    G4ThreeVector vertex2(xMax, yMax, zMin);
    G4ThreeVector vertex3(xMin, yMax, zMin);
    G4ThreeVector vertex4(xMin, yMin, zMax);
    G4ThreeVector vertex5(xMax, yMin, zMax);
    G4ThreeVector vertex6(xMax, yMax, zMax);
    G4ThreeVector vertex7(xMin, yMax, zMax);

    vertices->push_back(pTransform.TransformPoint(vertex0));
    vertices->push_back(pTransform.TransformPoint(vertex1));
    vertices->push_back(pTransform.TransformPoint(vertex2));
    vertices->push_back(pTransform.TransformPoint(vertex3));
    vertices->push_back(pTransform.TransformPoint(vertex4));
    vertices->push_back(pTransform.TransformPoint(vertex5));
    vertices->push_back(pTransform.TransformPoint(vertex6));
    vertices->push_back(pTransform.TransformPoint(vertex7));
  }
  else
  {
    G4Exception("G4VUSolid::CreateRotatedVertices()", "FatalError",
                FatalException, "Out of memory - Cannot allocate vertices!");
  }
  return vertices;
}

G4Polyhedron* G4USolid::CreatePolyhedron() const
{
  // Must be implemented in concrete wrappers...

  std::ostringstream message;
  message << "Visualization not supported for USolid shape "
          << GetEntityType() << "... Sorry!" << G4endl;
  G4Exception("G4USolid::CreatePolyhedron()", "GeomSolids0003",
              FatalException, message);
  return 0;
}

G4Polyhedron* G4USolid::GetPolyhedron() const
{
  if (!fPolyhedron ||
      fRebuildPolyhedron ||
      fPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fPolyhedron->GetNumberOfRotationSteps())
  {
    G4AutoLock l(&polyhedronMutex);
    delete fPolyhedron;
    fPolyhedron = CreatePolyhedron();
    fRebuildPolyhedron = false;
    l.unlock();
  }
  return fPolyhedron;
}

G4VisExtent G4USolid:: GetExtent() const
{
  G4VisExtent extent;
  G4VoxelLimits voxelLimits;  // Defaults to "infinite" limits.
  G4AffineTransform affineTransform;
  G4double vmin, vmax;
  CalculateExtent(kXAxis, voxelLimits, affineTransform, vmin, vmax);
  extent.SetXmin(vmin);
  extent.SetXmax(vmax);
  CalculateExtent(kYAxis, voxelLimits, affineTransform, vmin, vmax);
  extent.SetYmin(vmin);
  extent.SetYmax(vmax);
  CalculateExtent(kZAxis, voxelLimits, affineTransform, vmin, vmax);
  extent.SetZmin(vmin);
  extent.SetZmax(vmax);
  return extent;
}

#endif  // G4GEOM_USE_USOLIDS
