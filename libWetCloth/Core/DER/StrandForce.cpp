//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "StrandForce.h"

#include "../ThreadUtils.h"
#include "Dependencies/BendingProducts.h"
#include "Forces/BendingForce.h"
#include "Forces/ForceAccumulator.h"
#include "Forces/StretchingForce.h"
#include "Forces/TwistingForce.h"
#include "Forces/ViscousOrNotViscous.h"

// To match with rest of FilmFlow framework sign convention
//  (we compute Forces and Force Jacobians, FilmFlow expects Energy gradients
//  and Hessians)
#define FORCE_SIGN -1.0
#define HESS_SIGN -1.0

#define STRETCH
#define TWIST
#define BEND

StrandForce::StrandForce(const std::shared_ptr<TwoDScene>& scene,
                         const std::vector<int>& consecutiveStrandVertices,
                         const int& parameterIndex, int globalIndex)
    : m_verts(consecutiveStrandVertices),
      m_globalIndex(globalIndex),
      m_strandParams(NULL),
      m_scene(scene),
      m_strandEnergyUpdate(0.),
      m_strandForceUpdate(getNumVertices() * 4 - 1),
      m_strandHessianUpdate(),
      m_strandState(NULL),
      m_startState(NULL) {
  m_strandParams = m_scene->getElasticParameters(parameterIndex);

  VecX initDoFs(getNumVertices() * 4 - 1);
  for (int i = 0; i < getNumVertices(); ++i) {
    if (m_scene->isTip(m_verts[i]))
      initDoFs.segment<3>(i * 4) =
          m_scene->getX().segment<3>(m_scene->getDof(m_verts[i]));
    else
      initDoFs.segment<4>(i * 4) =
          m_scene->getX().segment<4>(m_scene->getDof(m_verts[i]));
  }
  m_strandState =
      new StrandState(initDoFs, m_strandParams->getBendingMatrixBase());
  m_startState = new StartState(initDoFs);

  m_packing_fraction.resize(m_verts.size());
  m_packing_fraction.setOnes();

  m_v_plus.resize(m_verts.size() * 4);
  m_v_plus.setZero();

  m_stretching_multipliers.resize(m_verts.size());
  m_bending_multipliers.resize(m_verts.size() * 4);
  m_twisting_multipliers.resize(m_verts.size());

  m_viscous_stretching_multipliers.resize(m_verts.size());
  m_viscous_bending_multipliers.resize(m_verts.size() * 4);
  m_viscous_twisting_multipliers.resize(m_verts.size());

  m_stretching_multipliers.setZero();
  m_bending_multipliers.setZero();
  m_twisting_multipliers.setZero();

  m_viscous_stretching_multipliers.setZero();
  m_viscous_bending_multipliers.setZero();
  m_viscous_twisting_multipliers.setZero();

  resizeInternals();
  freezeRestShape(
      0, getNumEdges());  // for now the rest shape is the shape in which the
                          // strand is created, unless this is called later on.

  // wont need to update first step's DoFs,
  // therefore must initialize the stored quantities as well
  recomputeGlobal();
}

StrandForce::~StrandForce() {}

StrandState::StrandState(const VecX& initDoFs,
                         BendingMatrixBase& bendingMatrixBase)
    : m_dofs(initDoFs),
      m_edges(m_dofs),
      m_lengths(m_edges),
      m_tangents(m_edges, m_lengths),
      m_referenceFrames1(m_tangents),
      m_referenceFrames2(m_tangents, m_referenceFrames1),
      m_referenceTwists(m_tangents, m_referenceFrames1),
      m_twists(m_referenceTwists, m_dofs),
      m_curvatureBinormals(m_tangents),
      m_trigThetas(m_dofs),
      m_materialFrames1(m_trigThetas, m_referenceFrames1, m_referenceFrames2),
      m_materialFrames2(m_trigThetas, m_referenceFrames1, m_referenceFrames2),
      m_kappas(m_curvatureBinormals, m_materialFrames1, m_materialFrames2),
      m_gradKappas(m_lengths, m_tangents, m_curvatureBinormals,
                   m_materialFrames1, m_materialFrames2, m_kappas),
      m_gradTwists(m_lengths, m_curvatureBinormals),
      m_gradTwistsSquared(m_gradTwists),
      m_hessKappas(m_lengths, m_tangents, m_curvatureBinormals,
                   m_materialFrames1, m_materialFrames2, m_kappas),
      m_hessTwists(m_tangents, m_lengths, m_curvatureBinormals),
      m_bendingProducts(bendingMatrixBase, m_gradKappas) {}

StartState::StartState(const VecX& initDoFs)
    : m_dofs(initDoFs),
      m_edges(m_dofs),
      m_lengths(m_edges),
      m_tangents(m_edges, m_lengths),
      m_referenceFrames1(m_tangents),
      m_referenceFrames2(m_tangents, m_referenceFrames1),
      m_referenceTwists(m_tangents, m_referenceFrames1),
      m_twists(m_referenceTwists, m_dofs),
      m_curvatureBinormals(m_tangents),
      m_trigThetas(m_dofs),
      m_materialFrames1(m_trigThetas, m_referenceFrames1, m_referenceFrames2),
      m_materialFrames2(m_trigThetas, m_referenceFrames1, m_referenceFrames2),
      m_kappas(m_curvatureBinormals, m_materialFrames1, m_materialFrames2) {}

void StrandForce::updateStartState() {
  const VectorXs& x = m_scene->getX();
  const VectorXs& psi = m_scene->getVolumeFraction();

  VecX currentStrandDoFs(getNumVertices() * 4);
  currentStrandDoFs.setZero();
  for (int i = 0; i < getNumVertices(); ++i) {
    if (m_scene->isTip(m_verts[i]))
      currentStrandDoFs.segment<3>(i * 4) =
          m_scene->getX().segment<3>(m_scene->getDof(m_verts[i]));
    else
      currentStrandDoFs.segment<4>(i * 4) =
          m_scene->getX().segment<4>(m_scene->getDof(m_verts[i]));
  }
  m_startState->m_dofs.set(currentStrandDoFs);

  const int num_verts = getNumVertices();
  for (int i = 0; i < num_verts; ++i) {
    m_packing_fraction[i] =
        pow(psi[m_verts[i]], m_scene->getSimInfo().lambda);
  }
}

void StrandForce::updateStrandState() {
  VecX initDoFs(getNumVertices() * 4);
  initDoFs.setZero();
  for (int i = 0; i < getNumVertices(); ++i) {
    if (m_scene->isTip(m_verts[i]))
      initDoFs.segment<3>(i * 4) =
          m_scene->getX().segment<3>(m_scene->getDof(m_verts[i]));
    else
      initDoFs.segment<4>(i * 4) =
          m_scene->getX().segment<4>(m_scene->getDof(m_verts[i]));
  }
  //    m_strandState->m_dofs.sets
  m_strandState->m_dofs.set(initDoFs);
}

void StrandForce::updateRestShape(const VecX& dof_restshape, scalar damping) {
  StartState restshape_state(dof_restshape);

  int nedges = getNumEdges();
  for (IndexType vtx = 0; vtx < nedges; ++vtx) {  // Fix rest lengths
    m_restLengths[vtx] = (1. - damping) * restshape_state.m_lengths[vtx] +
                         damping * m_restLengths[vtx];
  }
  updateEverythingThatDependsOnRestLengths();

  for (IndexType vtx = 0; vtx < nedges; ++vtx) {
    m_restKappas[vtx] = (1. - damping) * restshape_state.m_kappas[vtx] +
                        damping * m_restKappas[vtx];
    m_restTwists[vtx] = (1. - damping) * restshape_state.m_twists[vtx] +
                        damping * m_restTwists[vtx];
  }
}

void StrandForce::resizeInternals() {  // To be called on creation
  m_restLengths.resize(getNumEdges());
  m_restKappas.resize(getNumEdges());
  m_restTwists.resize(getNumEdges());
  m_vertexMasses.resize(getNumVertices());
  m_VoronoiLengths.resize(getNumVertices());
  m_invVoronoiLengths.resize(getNumVertices());
}

void StrandForce::freezeRestShape(
    unsigned begin, unsigned end,
    scalar damping) {  // Take the current configuration as rest shape

  for (IndexType vtx = begin; vtx < end; ++vtx) {  // Fix rest lengths
    m_restLengths[vtx] = (1. - damping) * m_strandState->m_lengths[vtx] +
                         damping * m_restLengths[vtx];
  }
  updateEverythingThatDependsOnRestLengths();

  for (IndexType vtx = begin; vtx < end; ++vtx) {
    m_restKappas[vtx] = (1. - damping) * m_strandState->m_kappas[vtx] +
                        damping * m_restKappas[vtx];
    m_restTwists[vtx] = (1. - damping) * m_strandState->m_twists[vtx] +
                        damping * m_restTwists[vtx];
  }
}

void StrandForce::updateEverythingThatDependsOnRestLengths() {
  // Total rest length
  m_totalRestLength = 0.0;
  for (IndexType vtx = 0; vtx < getNumEdges(); ++vtx) {
    m_totalRestLength += m_restLengths[vtx];
  }

  // Compute Voronoi lengths
  m_VoronoiLengths[0] = 0.5 * m_restLengths[0];
  for (IndexType vtx = 1; vtx < getNumEdges(); ++vtx) {
    m_VoronoiLengths[vtx] = 0.5 * (m_restLengths[vtx - 1] + m_restLengths[vtx]);
  }
  m_VoronoiLengths[getNumEdges()] = 0.5 * m_restLengths[getNumVertices() - 2];

  // Compute masses and inverse of Voronoi lengths
  for (IndexType vtx = 0; vtx < getNumVertices(); ++vtx) {
    m_vertexMasses[vtx] = m_strandParams->m_density * m_VoronoiLengths[vtx] *
                          M_PI * m_strandParams->getRadiusA(vtx) *
                          m_strandParams->getRadiusB(vtx);
    m_invVoronoiLengths[vtx] = 1.0 / m_VoronoiLengths[vtx];
  }
}

///////////////////////////////////////////// Force functions
/////////////////////////////////////////////////////////////

void StrandForce::clearStored() {
  m_strandEnergyUpdate = 0.;
  m_strandForceUpdate.setZero();
  m_strandHessianUpdate.clear();
  m_strandAngularHessianUpdate.clear();
}

void StrandForce::recomputeGlobal() {
  clearStored();
  accumulateQuantity(m_strandEnergyUpdate);
  accumulateQuantity(m_strandForceUpdate);
  accumulateHessian(m_strandHessianUpdate, m_strandAngularHessianUpdate);

  for (int i = 0; i < scripted_twist_nodes.size(); ++i)
  {
    ForceAccumulator<TwistingForce<NonViscous> >::accumulateScripted(
      m_strandForceUpdate, (IndexType)scripted_twist_nodes[i], (scalar)scripted_twist_forces[i]);
  }
  scripted_twist_nodes.clear();
  scripted_twist_forces.clear();

  // Free some memory
  m_strandState->m_hessTwists.free();
  m_strandState->m_hessKappas.free();
}

template <typename AccumulatedT>
void StrandForce::accumulateQuantity(AccumulatedT& accumulated) {
  ForceAccumulator<StretchingForce<NonViscous> >::accumulate(accumulated,
                                                             *this);
  ForceAccumulator<TwistingForce<NonViscous> >::accumulate(accumulated, *this);
  ForceAccumulator<BendingForce<NonViscous> >::accumulate(accumulated, *this);

  if (m_strandParams->m_accumulateWithViscous) {
    if (!m_strandParams->m_accumulateViscousOnlyForBendingModes) {
      ForceAccumulator<StretchingForce<Viscous> >::accumulate(accumulated,
                                                              *this);
    }
    ForceAccumulator<TwistingForce<Viscous> >::accumulate(accumulated, *this);
    ForceAccumulator<BendingForce<Viscous> >::accumulate(accumulated, *this);
  }
}

void StrandForce::accumulateHessian(TripletXs& accumulated,
                                    TripletXs& accumulated_twist) {
  ForceAccumulator<StretchingForce<NonViscous> >::accumulate(
      accumulated, accumulated_twist, *this);
  ForceAccumulator<TwistingForce<NonViscous> >::accumulate(
      accumulated, accumulated_twist, *this);
  ForceAccumulator<BendingForce<NonViscous> >::accumulate(
      accumulated, accumulated_twist, *this);

  if (m_strandParams->m_accumulateWithViscous) {
    if (!m_strandParams->m_accumulateViscousOnlyForBendingModes) {
      ForceAccumulator<StretchingForce<Viscous> >::accumulate(
          accumulated, accumulated_twist, *this);
    }
    ForceAccumulator<TwistingForce<Viscous> >::accumulate(
        accumulated, accumulated_twist, *this);
    ForceAccumulator<BendingForce<Viscous> >::accumulate(
        accumulated, accumulated_twist, *this);
  }
}

Force* StrandForce::createNewCopy() { return new StrandForce(*this); }

void StrandForce::preCompute() {
  /* nothing to do here, updateStartDoFs called separately and otherwise need to
   * update every time we compute (in case nonlinear) */
  updateStrandState();
  recomputeGlobal();
}

int StrandForce::numConstraintNonViscous() {  // Spring = numEdges  //Twist =
                                              // NumVertices - 2  //Bending =
                                              // 4*(NumVertices - 2)
  int numConstraintNonViscous = 0;
#ifdef STRETCH
  numConstraintNonViscous += getNumEdges();  // Stretch
#endif

#ifdef TWIST
  numConstraintNonViscous += (getNumVertices() - 2);  // Twist
#endif

#ifdef BEND
  numConstraintNonViscous += 4 * (getNumVertices() - 2);  // Bend
#endif
  return numConstraintNonViscous;
}

int StrandForce::numConstraintViscous() {
  int numConstraintViscous = 0;
  if (m_strandParams->m_accumulateWithViscous) {
    if (!m_strandParams->m_accumulateViscousOnlyForBendingModes) {
#ifdef STRETCH
      numConstraintViscous += getNumEdges();  // Stretch
#endif
    }
#ifdef TWIST
    numConstraintViscous += (getNumVertices() - 2);  // Twist
#endif
#ifdef BEND
    numConstraintViscous += 4 * (getNumVertices() - 2);  // Bend
#endif
  }
  return numConstraintViscous;
}

const char* StrandForce::name() { return "Strand Material Forces"; }

void StrandForce::addEnergyToTotal(const VectorXs& x, const VectorXs& v,
                                   const VectorXs& m, const VectorXs& psi,
                                   const scalar& lambda, scalar& E) {
  // TODO
  E += m_strandEnergyUpdate;
}

void StrandForce::addGradEToTotal(const VectorXs& x, const VectorXs& v,
                                  const VectorXs& m, const VectorXs& psi,
                                  const scalar& lambda, VectorXs& gradE) {
  const int num_verts = m_verts.size();

  threadutils::for_each(0, num_verts, [&](int i) {
    if (i != num_verts - 1)
      gradE.segment<4>(4 * m_verts[i]) -= m_strandForceUpdate.segment<4>(i * 4);
    else
      gradE.segment<3>(4 * m_verts[i]) -= m_strandForceUpdate.segment<3>(i * 4);
  });
}

void StrandForce::addHessXToTotal(const VectorXs& x, const VectorXs& v,
                                  const VectorXs& m, const VectorXs& psi,
                                  const scalar& lambda, TripletXs& hessE,
                                  int hessE_index, const scalar& dt) {
  const int num_hess = numHessX();

  threadutils::for_each(0, num_hess, [&](int i) {
    const Triplets& data = m_strandHessianUpdate[i];
    int col_vert = data.col() / 4;
    int col_r = data.col() - col_vert * 4;
    int row_vert = data.row() / 4;
    int row_r = data.row() - row_vert * 4;
    hessE[hessE_index + i] =
        Triplets(4 * m_verts[row_vert] + row_r, 4 * m_verts[col_vert] + col_r,
                 -data.value());
  });
}

void StrandForce::addAngularHessXToTotal(const VectorXs& x, const VectorXs& v,
                                         const VectorXs& m, const VectorXs& psi,
                                         const scalar& lambda, TripletXs& hessE,
                                         int hessE_index, const scalar& dt) {
  const int num_hess = numAngularHessX();

  threadutils::for_each(0, num_hess, [&](int i) {
    const Triplets& data = m_strandAngularHessianUpdate[i];
    int col_vert = data.col() / 4;
    int row_vert = data.row() / 4;
    hessE[hessE_index + i] =
        Triplets(m_verts[row_vert], m_verts[col_vert], -data.value());
  });
}

void StrandForce::updateMultipliers(const VectorXs& x, const VectorXs& vplus,
                                    const VectorXs& m, const VectorXs& psi,
                                    const scalar& lambda, const scalar& dt) {
  if (m_strandParams->m_useApproxJacobian ||
      !m_strandParams->m_useTournierJacobian) {
    return;
  }

  const int num_verts = getNumVertices();
  for (int i = 0; i < num_verts; ++i) {
    m_v_plus.segment<4>(i * 4) = vplus.segment<4>(m_verts[i] * 4);
  }

  m_stretching_multipliers.setZero();
  m_twisting_multipliers.setZero();
  m_bending_multipliers.setZero();

  ForceAccumulator<StretchingForce<NonViscous> >::accumulateMultipliers(
      m_stretching_multipliers, *this, dt);
  ForceAccumulator<TwistingForce<NonViscous> >::accumulateMultipliers(
      m_twisting_multipliers, *this, dt);
  ForceAccumulator<BendingForce<NonViscous> >::accumulateMultipliers(
      m_bending_multipliers, *this, dt);

  if (m_strandParams->m_accumulateWithViscous) {
    if (!m_strandParams->m_accumulateViscousOnlyForBendingModes) {
      m_viscous_stretching_multipliers.setZero();
      ForceAccumulator<StretchingForce<Viscous> >::accumulateMultipliers(
          m_viscous_stretching_multipliers, *this, dt);
    }
    m_viscous_twisting_multipliers.setZero();
    m_viscous_bending_multipliers.setZero();

    ForceAccumulator<TwistingForce<Viscous> >::accumulateMultipliers(
        m_viscous_twisting_multipliers, *this, dt);
    ForceAccumulator<BendingForce<Viscous> >::accumulateMultipliers(
        m_viscous_bending_multipliers, *this, dt);
  }
}

int StrandForce::numHessX() { return m_strandHessianUpdate.size(); }

int StrandForce::numAngularHessX() {
  return m_strandAngularHessianUpdate.size();
}

int StrandForce::flag() const { return 1; }

void StrandForce::AddScriptedForce(int node, scalar t)
{
  scripted_twist_nodes.push_back(node);
  scripted_twist_forces.push_back(t);
}
