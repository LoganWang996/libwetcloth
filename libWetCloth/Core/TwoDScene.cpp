//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TwoDScene.h"

#include <igl/point_simplex_squared_distance.h>
#include <igl/ray_mesh_intersect.h>

#include <iostream>
#include <numeric>
#include <stack>

#include "AttachForce.h"
#include "DER/StrandForce.h"
#include "MathUtilities.h"
#include "ThreadUtils.h"
#include "SpherePattern.h"

/*!
 * Outputting parameters for debugging
 */
std::ostream& operator<<(std::ostream& os, const SimInfo& info) {
  os << "rest volume fraction: " << info.rest_volume_fraction << std::endl;
  os << "lambda: " << info.lambda << std::endl;
  os << "elasto flip stretching coeff: " << info.elasto_flip_coeff << std::endl;
  os << "elasto flip-asym coeff: " << info.elasto_flip_asym_coeff << std::endl;
  os << "elasto advection coeff: " << info.elasto_advect_coeff << std::endl;
  os << "solve solid: " << info.solve_solid << std::endl;
  os << "use twist: " << info.use_twist << std::endl;
  os << "use amgpcg solid: " << info.use_amgpcg_solid << std::endl;
  return os;
}

/*!
 * init the Scene structure
 */
TwoDScene::TwoDScene()
    : m_x(),
      m_v(),
      m_m(),
      m_fixed(),
      step_count(0),
      m_edges(),
      m_num_colors(1),
      m_forces() {
  sphere_pattern::generateSpherePattern(m_sphere_pattern);
}

TwoDScene::~TwoDScene() {}

int TwoDScene::getNumBuckets() const { return m_particle_buckets.size(); }

int TwoDScene::getNumParticles() const { return m_x.size() / 4; }

int TwoDScene::getNumEdges() const { return m_edges.rows(); }

int TwoDScene::getNumSurfels() const { return m_surfels.size(); }

int TwoDScene::getNumFaces() const { return m_faces.rows(); }

int TwoDScene::getNumGausses() const { return m_x_gauss.size() / 4; }

int TwoDScene::getNumElasticParameters() const {
  return m_strandParameters.size();
}

const VectorXs& TwoDScene::getX() const { return m_x; }

VectorXs& TwoDScene::getX() { return m_x; }

const VectorXs& TwoDScene::getV() const { return m_v; }

VectorXs& TwoDScene::getV() { return m_v; }

const VectorXs& TwoDScene::getFluidV() const { return m_fluid_v; }

VectorXs& TwoDScene::getFluidV() { return m_fluid_v; }

const VectorXs& TwoDScene::getM() const { return m_m; }

VectorXs& TwoDScene::getM() { return m_m; }

const VectorXs& TwoDScene::getFluidM() const { return m_fluid_m; }

VectorXs& TwoDScene::getFluidM() { return m_fluid_m; }

const VectorXs& TwoDScene::getFluidVol() const { return m_fluid_vol; }

VectorXs& TwoDScene::getFluidVol() { return m_fluid_vol; }

const VectorXs& TwoDScene::getVol() const { return m_vol; }

VectorXs& TwoDScene::getVol() { return m_vol; }

const VectorXs& TwoDScene::getRadius() const { return m_radius; }

VectorXs& TwoDScene::getRadius() { return m_radius; }

const MatrixXs& TwoDScene::getGaussFe() const { return m_Fe_gauss; }

MatrixXs& TwoDScene::getGaussFe() { return m_Fe_gauss; }

const VectorXs& TwoDScene::getGaussX() const { return m_x_gauss; }

VectorXs& TwoDScene::getGaussX() { return m_x_gauss; }

const VectorXs& TwoDScene::getGaussV() const { return m_v_gauss; }

VectorXs& TwoDScene::getGaussV() { return m_v_gauss; }

const VectorXs& TwoDScene::getGaussM() const { return m_m_gauss; }

VectorXs& TwoDScene::getGaussM() { return m_m_gauss; }

const VectorXs& TwoDScene::getGaussVol() const { return m_vol_gauss; }

VectorXs& TwoDScene::getGaussVol() { return m_vol_gauss; }

const MatrixXs& TwoDScene::getGaussd() const { return m_d_gauss; }

MatrixXs& TwoDScene::getGaussd() { return m_d_gauss; }

const MatrixXs& TwoDScene::getGaussDinv() const { return m_D_inv_gauss; }

MatrixXs& TwoDScene::getGaussDinv() { return m_D_inv_gauss; }

const MatrixXs& TwoDScene::getGaussD() const { return m_D_gauss; }

MatrixXs& TwoDScene::getGaussD() { return m_D_gauss; }

/*!
 * swap particle pos in buffer, used for particle deletion
 */
void TwoDScene::swapParticles(int i, int j) {
  mathutils::swap<scalar, 4>(m_x, i, j);
  mathutils::swap<scalar, 4>(m_rest_x, i, j);
  mathutils::swap<scalar, 4>(m_v, i, j);
  mathutils::swap<scalar, 4>(m_saved_v, i, j);
  mathutils::swap<scalar, 4>(m_dv, i, j);
  mathutils::swap<scalar, 4>(m_fluid_v, i, j);
  mathutils::swap<scalar, 4>(m_m, i, j);
  mathutils::swap<scalar, 4>(m_fluid_m, i, j);
  mathutils::swap<scalar, 3>(m_orientation, i, j);
  mathutils::swap<scalar, 2>(m_radius, i, j);
  std::swap(m_vol(i), m_vol(j));
  std::swap(m_rest_vol(i), m_rest_vol(j));
  std::swap(m_fluid_vol(i), m_fluid_vol(j));
  std::swap(m_shape_factor(i), m_shape_factor(j));
  std::swap(m_fixed[i], m_fixed[j]);
  std::swap(m_twist[i], m_twist[j]);
  std::swap(m_particle_to_edge[i], m_particle_to_edge[j]);
  std::swap(m_particle_to_face[i], m_particle_to_face[j]);
  std::swap(m_particle_to_surfel[i], m_particle_to_surfel[j]);
  std::swap(m_particle_rest_length[i], m_particle_rest_length[j]);
  std::swap(m_particle_rest_area[i], m_particle_rest_area[j]);
  std::swap(m_particle_group[i], m_particle_group[j]);
  std::swap(m_volume_fraction[i], m_volume_fraction[j]);
  std::swap(m_rest_volume_fraction[i], m_rest_volume_fraction[j]);
  std::swap(m_inside[i], m_inside[j]);
  std::swap(m_classifier[i], m_classifier[j]);

  std::swap(m_is_strand_tip[i], m_is_strand_tip[j]);

  mathutils::swap<scalar, 3>(m_B, i, j);
  mathutils::swap<scalar, 3>(m_fB, i, j);
}

const Matrix27x4s& TwoDScene::getParticleWeights(int pidx) const {
  return m_particle_weights[pidx];
}

const std::vector<unsigned char>& TwoDScene::getFixed() const {
  return m_fixed;
}

const std::vector<bool>& TwoDScene::getTwist() const { return m_twist; }

const std::vector<Matrix27x4s>& TwoDScene::getParticleWeights() const {
  return m_particle_weights;
}

Matrix27x4s& TwoDScene::getParticleWeights(int pidx) {
  return m_particle_weights[pidx];
}

const Matrix27x3s& TwoDScene::getGaussWeights(int pidx) const {
  return m_gauss_weights[pidx];
}

Matrix27x3s& TwoDScene::getGaussWeights(int pidx) {
  return m_gauss_weights[pidx];
}

const Matrix27x2i& TwoDScene::getParticleNodesSolidPhi(int pidx) const {
  return m_particle_nodes_solid_phi[pidx];
}

Matrix27x2i& TwoDScene::getParticleNodesSolidPhi(int pidx) {
  return m_particle_nodes_solid_phi[pidx];
}

const Matrix27x2i& TwoDScene::getParticleNodesX(int pidx) const {
  return m_particle_nodes_x[pidx];
}

Matrix27x2i& TwoDScene::getParticleNodesX(int pidx) {
  return m_particle_nodes_x[pidx];
}

const Matrix27x2i& TwoDScene::getParticleNodesY(int pidx) const {
  return m_particle_nodes_y[pidx];
}

Matrix27x2i& TwoDScene::getParticleNodesY(int pidx) {
  return m_particle_nodes_y[pidx];
}

const Matrix27x2i& TwoDScene::getParticleNodesZ(int pidx) const {
  return m_particle_nodes_z[pidx];
}

Matrix27x2i& TwoDScene::getParticleNodesZ(int pidx) {
  return m_particle_nodes_z[pidx];
}

const Matrix27x2i& TwoDScene::getGaussNodesX(int pidx) const {
  return m_gauss_nodes_x[pidx];
}

Matrix27x2i& TwoDScene::getGaussNodesX(int pidx) {
  return m_gauss_nodes_x[pidx];
}

const Matrix27x2i& TwoDScene::getGaussNodesY(int pidx) const {
  return m_gauss_nodes_y[pidx];
}

Matrix27x2i& TwoDScene::getGaussNodesY(int pidx) {
  return m_gauss_nodes_y[pidx];
}

const Matrix27x2i& TwoDScene::getGaussNodesZ(int pidx) const {
  return m_gauss_nodes_z[pidx];
}

Matrix27x2i& TwoDScene::getGaussNodesZ(int pidx) {
  return m_gauss_nodes_z[pidx];
}

int TwoDScene::getDefaultNumNodes() const { return m_num_nodes; }

int TwoDScene::getNumNodes(int bucket_idx) const {
  if (m_bucket_activated[bucket_idx])
    return m_num_nodes * m_num_nodes * m_num_nodes;
  else
    return 0;
}

const std::vector<VectorXuc>& TwoDScene::getNodeStateX() const {
  return m_node_state_u;
}

std::vector<VectorXuc>& TwoDScene::getNodeStateX() { return m_node_state_u; }

const std::vector<VectorXuc>& TwoDScene::getNodeStateY() const {
  return m_node_state_v;
}

std::vector<VectorXuc>& TwoDScene::getNodeStateY() { return m_node_state_v; }

const std::vector<VectorXuc>& TwoDScene::getNodeStateZ() const {
  return m_node_state_w;
}

std::vector<VectorXuc>& TwoDScene::getNodeStateZ() { return m_node_state_w; }

const std::vector<VectorXs>& TwoDScene::getNodeSolidPhi() const {
  return m_node_solid_phi;
}

std::vector<VectorXs>& TwoDScene::getNodeSolidPhi() { return m_node_solid_phi; }

const VectorXs& TwoDScene::getNodePos(int bucket_idx) const {
  return m_node_pos[bucket_idx];
}

VectorXs& TwoDScene::getNodePos(int bucket_idx) {
  return m_node_pos[bucket_idx];
}

const std::vector<VectorXs>& TwoDScene::getNodePos() const {
  return m_node_pos;
}

std::vector<VectorXs>& TwoDScene::getNodePos() { return m_node_pos; }

const std::vector<VectorXs>& TwoDScene::getNodeCellSolidPhi() const {
  return m_node_cell_solid_phi;
}

std::vector<VectorXs>& TwoDScene::getNodeCellSolidPhi() {
  return m_node_cell_solid_phi;
}

const std::vector<VectorXs>& TwoDScene::getNodePressure() const {
  return m_node_pressure;
}

std::vector<VectorXs>& TwoDScene::getNodePressure() { return m_node_pressure; }

const std::vector<VectorXs>& TwoDScene::getNodeVelocityX() const {
  return m_node_vel_x;
}

std::vector<VectorXs>& TwoDScene::getNodeVelocityX() { return m_node_vel_x; }

const std::vector<VectorXs>& TwoDScene::getNodeVelocityY() const {
  return m_node_vel_y;
}

std::vector<VectorXs>& TwoDScene::getNodeVelocityY() { return m_node_vel_y; }

const std::vector<VectorXs>& TwoDScene::getNodeVelocityZ() const {
  return m_node_vel_z;
}

std::vector<VectorXs>& TwoDScene::getNodeVelocityZ() { return m_node_vel_z; }

const std::vector<VectorXs>& TwoDScene::getNodeFluidVelocityX() const {
  return m_node_vel_fluid_x;
}

std::vector<VectorXs>& TwoDScene::getNodeFluidVelocityX() {
  return m_node_vel_fluid_x;
}

const std::vector<VectorXs>& TwoDScene::getNodeFluidVelocityY() const {
  return m_node_vel_fluid_y;
}

std::vector<VectorXs>& TwoDScene::getNodeFluidVelocityY() {
  return m_node_vel_fluid_y;
}

const std::vector<VectorXs>& TwoDScene::getNodeFluidVelocityZ() const {
  return m_node_vel_fluid_z;
}

std::vector<VectorXs>& TwoDScene::getNodeFluidVelocityZ() {
  return m_node_vel_fluid_z;
}

const std::vector<unsigned char>& TwoDScene::getBucketActivated() const {
  return m_bucket_activated;
}

const std::vector<VectorXs>& TwoDScene::getNodeMassX() const {
  return m_node_mass_x;
}

std::vector<VectorXs>& TwoDScene::getNodeMassX() { return m_node_mass_x; }

const std::vector<VectorXs>& TwoDScene::getNodeMassY() const {
  return m_node_mass_y;
}

std::vector<VectorXs>& TwoDScene::getNodeMassY() { return m_node_mass_y; }

const std::vector<VectorXs>& TwoDScene::getNodeMassZ() const {
  return m_node_mass_z;
}

std::vector<VectorXs>& TwoDScene::getNodeMassZ() { return m_node_mass_z; }

const std::vector<VectorXs>& TwoDScene::getNodeVolX() const {
  return m_node_vol_x;
}

std::vector<VectorXs>& TwoDScene::getNodeVolX() { return m_node_vol_x; }

const std::vector<VectorXs>& TwoDScene::getNodeVolY() const {
  return m_node_vol_y;
}

std::vector<VectorXs>& TwoDScene::getNodeVolY() { return m_node_vol_y; }

const std::vector<VectorXs>& TwoDScene::getNodeVolZ() const {
  return m_node_vol_z;
}

std::vector<VectorXs>& TwoDScene::getNodeVolZ() { return m_node_vol_z; }

/*!
 * mark the inside/outside of a levelset
 */
void TwoDScene::markInsideOut() {
  const int num_parts = getNumParticles();
  threadutils::for_each(0, num_parts, [this](int pidx) {
    bool has_compressed = false;
    bool has_uncompressed = false;

    const auto& indices_x = m_particle_nodes_x[pidx];
    const int num_rows_x = indices_x.rows();
    for (int i = 0; i < num_rows_x; ++i) {
      if (m_bucket_activated[indices_x(i, 0)]) {
        has_compressed = true;
      } else {
        has_uncompressed = true;
      }
    }

    const auto& indices_y = m_particle_nodes_y[pidx];
    const int num_rows_y = indices_y.rows();
    for (int i = 0; i < num_rows_y; ++i) {
      if (m_bucket_activated[indices_y(i, 0)]) {
        has_compressed = true;
      } else {
        has_uncompressed = true;
      }
    }

    const auto& indices_z = m_particle_nodes_z[pidx];
    const int num_rows_z = indices_z.rows();
    for (int i = 0; i < num_rows_z; ++i) {
      if (m_bucket_activated[indices_z(i, 0)]) {
        has_compressed = true;
      } else {
        has_uncompressed = true;
      }
    }

    if (has_compressed) {
      if (has_uncompressed)
        m_inside[pidx] = 1U;
      else
        m_inside[pidx] = 2U;
    } else {
      m_inside[pidx] = 0U;
    }
  });
}

const VectorXuc& TwoDScene::getOutsideInfo() const { return m_inside; }

const std::vector<VectorXs>& TwoDScene::getNodeFluidVolX() const {
  return m_node_vol_fluid_x;
}

std::vector<VectorXs>& TwoDScene::getNodeFluidVolX() {
  return m_node_vol_fluid_x;
}

const std::vector<VectorXs>& TwoDScene::getNodeFluidVolY() const {
  return m_node_vol_fluid_y;
}

std::vector<VectorXs>& TwoDScene::getNodeFluidVolY() {
  return m_node_vol_fluid_y;
}

const std::vector<VectorXs>& TwoDScene::getNodeFluidVolZ() const {
  return m_node_vol_fluid_z;
}

std::vector<VectorXs>& TwoDScene::getNodeFluidVolZ() {
  return m_node_vol_fluid_z;
}

const std::vector<VectorXs>& TwoDScene::getNodeFluidMassX() const {
  return m_node_mass_fluid_x;
}

std::vector<VectorXs>& TwoDScene::getNodeFluidMassX() {
  return m_node_mass_fluid_x;
}

const std::vector<VectorXs>& TwoDScene::getNodeFluidMassY() const {
  return m_node_mass_fluid_y;
}

std::vector<VectorXs>& TwoDScene::getNodeFluidMassY() {
  return m_node_mass_fluid_y;
}

const std::vector<VectorXs>& TwoDScene::getNodeFluidMassZ() const {
  return m_node_mass_fluid_z;
}

std::vector<VectorXs>& TwoDScene::getNodeFluidMassZ() {
  return m_node_mass_fluid_z;
}

scalar TwoDScene::getFrictionAlpha(int pidx) const {
  return m_strandParameters[m_gauss_to_parameters[pidx]]->m_friction_alpha;
}

scalar TwoDScene::getFrictionBeta(int pidx) const {
  return m_strandParameters[m_gauss_to_parameters[pidx]]->m_friction_beta;
}

scalar TwoDScene::getGaussDensity(int pidx) const {
  return m_strandParameters[m_gauss_to_parameters[pidx]]->m_density;
}

scalar TwoDScene::getInitialVolumeFraction(int pidx) const {
  return m_sim_info.rest_volume_fraction;
}

/*!
 * get the equivalent radius for Elements
 */
scalar TwoDScene::getGaussRadius(int pidx, int dir) const {
  const int num_edges = getNumEdges();
  const int num_faces = getNumFaces();

  const scalar dx = getCellSize() * 0.125;

  if (pidx < num_edges) {
    const scalar r0 = m_radius[m_edges(pidx, 0) * 2 + dir];
    const scalar r1 = m_radius[m_edges(pidx, 1) * 2 + dir];
    return sqrt((r0 * r0 + r1 * r1) * 0.5);
  } else if (pidx < num_edges + num_faces) {
    const scalar r0 = m_radius[m_faces(pidx - num_edges, 0) * 2 + dir];
    const scalar r1 = m_radius[m_faces(pidx - num_edges, 1) * 2 + dir];
    const scalar r2 = m_radius[m_faces(pidx - num_edges, 2) * 2 + dir];
    return sqrt((r0 * r0 + r1 * r1 + r2 * r2) / 3.0);
  } else {
    return mathutils::defaultRadiusMultiplier() * dx * 8.0;
  }
}

scalar TwoDScene::getMu(int pidx) const {
  return m_strandParameters[m_gauss_to_parameters[pidx]]->m_shearModulus.get();
}

scalar TwoDScene::getLa(int pidx) const {
  double mu, la, E;
  mu = m_strandParameters[m_gauss_to_parameters[pidx]]->m_shearModulus.get();
  E = m_strandParameters[m_gauss_to_parameters[pidx]]->m_youngsModulus.get();
  la = mu * (E - 2 * mu) / (3 * mu - E);
  return la;
}

scalar TwoDScene::getAttachMultiplier(int pidx) const {
  return m_strandParameters[m_gauss_to_parameters[pidx]]->m_attachMultiplier;
}

scalar TwoDScene::getYoungModulus(int pidx) const {
  return m_strandParameters[m_gauss_to_parameters[pidx]]->m_youngsModulus.get();
}

scalar TwoDScene::getViscousModulus(int pidx) const {
  return m_strandParameters[m_gauss_to_parameters[pidx]]->m_viscosity;
}

scalar TwoDScene::getShearModulus(int pidx) const {
  return m_strandParameters[m_gauss_to_parameters[pidx]]->m_shearModulus.get();
}

scalar TwoDScene::getCollisionMultiplier(int pidx) const {
  return m_strandParameters[m_gauss_to_parameters[pidx]]->m_collisionMultiplier;
}

/*!
 * load forces fixing vertices to some specific position
 */
void TwoDScene::loadAttachForces() {
  const int num_part = getNumSoftElastoParticles();
  const int num_edges = getNumEdges();

  for (int i = 0; i < num_part; ++i) {
    if (m_particle_to_surfel[i] >= 0) continue;

    scalar K = 0.0;
    scalar mu = 0.0;
    scalar w = 0.0;

    for (int e : m_particle_to_edge[i]) {
      K += getYoungModulus(e) * getAttachMultiplier(e) * 0.5 * m_vol_gauss(e);
      mu += getShearModulus(e) * getAttachMultiplier(e) * 0.5 * m_vol_gauss(e);
      w += 0.5 * m_vol_gauss(e);
    }

    for (auto& p : m_particle_to_face[i]) {
      K += getYoungModulus(p.first + num_edges) *
           getAttachMultiplier(p.first + num_edges) * p.second *
           m_vol_gauss(p.first + num_edges);
      mu += getShearModulus(p.first + num_edges) *
            getAttachMultiplier(p.first + num_edges) * p.second *
            m_vol_gauss(p.first + num_edges);
      w += p.second * m_vol_gauss(p.first + num_edges);
    }

    if (w < 1e-20) continue;

    K /= w;
    mu /= w;

    scalar ks = 0.0;
    scalar kt = 0.0;
    scalar bs = 0.0;
    scalar bt = 0.0;
    if (isFixed(i) & 1) {
      ks = K * pow(m_rest_vol(i), 1. / 3.);
    }

    if (m_twist[i] && (isFixed(i) & 2)) {
      kt = mu * M_PI_4 * m_radius(i * 2) * m_radius(i * 2 + 1) *
           (m_radius(i * 2) * m_radius(i * 2) +
            m_radius(i * 2 + 1) * m_radius(i * 2 + 1)) /
           m_particle_rest_length(i);
    }

    if (ks > 0.0 || kt > 0.0) {
      std::shared_ptr<AttachForce> af =
          std::make_shared<AttachForce>(i, shared_from_this(), ks, kt, bs, bt);
      m_forces.push_back(af);
      m_attach_forces.push_back(af);
    }
  }
}

bool TwoDScene::isTip(int particle) const {
  assert(particle >= 0);
  assert(particle < getNumParticles());

  return m_is_strand_tip[particle];
}

bool TwoDScene::isSoft(int pidx) const {
  return m_particle_to_surfel[pidx] < 0;
}

/*!
 * resize the system
 */
void TwoDScene::resizeParticleSystem(int num_particles) {
  assert(num_particles >= 0);

  m_x.resize(4 * num_particles);
  m_rest_x.resize(4 * num_particles);
  m_v.resize(4 * num_particles);
  m_saved_v.resize(4 * num_particles);
  m_dv.resize(4 * num_particles);
  m_fluid_v.resize(4 * num_particles);
  m_m.resize(4 * num_particles);
  m_fluid_m.resize(4 * num_particles);
  m_vol.resize(num_particles);
  m_rest_vol.resize(num_particles);
  m_shape_factor.resize(num_particles);
  m_fluid_vol.resize(num_particles);
  m_radius.resize(2 * num_particles);
  m_fixed.resize(num_particles);
  m_twist.resize(num_particles);
  m_particle_to_edge.resize(num_particles);
  m_particle_to_surfel.resize(num_particles, -1);
  m_particle_to_face.resize(num_particles);
  m_particle_rest_length.resize(num_particles);
  m_particle_rest_area.resize(num_particles);
  m_particle_group.resize(num_particles);
  m_volume_fraction.resize(num_particles);
  m_rest_volume_fraction.resize(num_particles);
  m_div.resize(num_particles);
  m_inside.resize(num_particles);
  m_classifier.resize(num_particles, PC_NONE);
  m_orientation.resize(3 * num_particles);

  m_B.resize(num_particles * 3, 3);
  m_fB.resize(num_particles * 3, 3);

  m_particle_nodes_x.resize(num_particles);
  m_particle_nodes_y.resize(num_particles);
  m_particle_nodes_z.resize(num_particles);
  m_particle_nodes_p.resize(num_particles);
  m_particle_nodes_solid_phi.resize(num_particles);

  m_particle_weights.resize(num_particles);
  m_particle_weights_p.resize(num_particles);

  m_is_strand_tip.resize(num_particles);

  m_particle_rest_length.setZero();
  m_particle_rest_area.setZero();
  m_x.setZero();
  m_v.setZero();
  m_saved_v.setZero();
  m_dv.setZero();
  m_fluid_v.setZero();
  m_m.setZero();
  m_fluid_m.setZero();
  m_vol.setOnes();
  m_rest_vol.setOnes();
  m_fluid_vol.setZero();
  m_radius.setOnes();
  m_B.setZero();
  m_fB.setZero();
  m_volume_fraction.setZero();
  m_rest_volume_fraction.setZero();
  m_inside.setZero();
}

bool TwoDScene::isOutsideFluid(int pidx) const {
  return isFluid(pidx) && m_inside[pidx] == 0U;
}

/*!
 * resize the system in a conservative way
 */
void TwoDScene::conservativeResizeParticles(int num_particles) {
  m_x.conservativeResize(4 * num_particles);
  m_rest_x.conservativeResize(4 * num_particles);
  m_v.conservativeResize(4 * num_particles);
  m_dv.conservativeResize(4 * num_particles);
  m_saved_v.conservativeResize(4 * num_particles);
  m_fluid_v.conservativeResize(4 * num_particles);
  m_m.conservativeResize(4 * num_particles);
  m_fluid_m.conservativeResize(4 * num_particles);
  m_fluid_vol.conservativeResize(num_particles);
  m_vol.conservativeResize(num_particles);
  m_rest_vol.conservativeResize(num_particles);
  m_shape_factor.conservativeResize(num_particles);
  m_radius.conservativeResize(2 * num_particles);
  m_volume_fraction.conservativeResize(num_particles);
  m_rest_volume_fraction.conservativeResize(num_particles);
  m_fixed.resize(num_particles);
  m_twist.resize(num_particles);
  m_particle_to_edge.resize(num_particles);
  m_particle_to_face.resize(num_particles);
  m_particle_to_surfel.resize(num_particles);
  m_particle_rest_length.conservativeResize(num_particles);
  m_particle_rest_area.conservativeResize(num_particles);
  m_particle_group.resize(num_particles);
  m_inside.resize(num_particles);
  m_classifier.resize(num_particles);
  m_orientation.conservativeResize(3 * num_particles);

  m_B.conservativeResize(num_particles * 3, 3);
  m_fB.conservativeResize(num_particles * 3, 3);

  m_particle_nodes_x.resize(num_particles);
  m_particle_nodes_y.resize(num_particles);
  m_particle_nodes_z.resize(num_particles);
  m_particle_nodes_p.resize(num_particles);
  m_particle_nodes_solid_phi.resize(num_particles);

  m_particle_weights.resize(num_particles);
  m_particle_weights_p.resize(num_particles);

  m_is_strand_tip.resize(num_particles);
  m_div.resize(num_particles);
}

/*!
 * resize edges
 */
void TwoDScene::conservativeResizeEdges(int num_edges) {
  m_edges.conservativeResize(num_edges, 2);
  m_edge_rest_length.conservativeResize(num_edges);
  m_edge_inv_mapping.resize(num_edges);
  m_edge_to_parameters.resize(num_edges);
}

void TwoDScene::setEdgeToParameter(int idx, int params) {
  m_edge_to_parameters[idx] = params;
}

void TwoDScene::setFaceToParameter(int idx, int params) {
  m_face_to_parameters[idx] = params;
}

/*!
 * resize face elements
 */
void TwoDScene::conservativeResizeFaces(int num_faces) {
  m_faces.conservativeResize(num_faces, 3);
  m_face_weights.resize(num_faces);
  m_face_rest_area.conservativeResize(num_faces);
  m_face_inv_mapping.resize(num_faces);
  m_face_to_parameters.resize(num_faces);
}

/*!
 * calculate local divergence on particles
 */
void TwoDScene::updateParticleDiv() {
  const int num_particles = getNumParticles();

  // update connectivity
  const int num_edges = m_edges.rows();
  const int num_triangles = m_faces.rows();

  threadutils::for_each(0, num_particles, [&](int pidx) {
    const int num_ele =
        m_particle_to_edge[pidx].size() + m_particle_to_face[pidx].size();
    m_div[pidx].resize(num_ele * 3);
    m_div[pidx].setZero();
  });

  m_gauss_buckets.for_each_bucket_particles_colored([&](int gidx, int) {
    if (gidx < num_edges) {
      const int eidx = gidx;
      const Vector2i& im = m_edge_inv_mapping[eidx];
      const auto& e = m_edges.row(eidx);

      Vector3s ev = m_x.segment<3>(e(1) * 4) - m_x.segment<3>(e(0) * 4);
      const scalar evl = ev.norm();
      if (evl > 1e-20) {
        ev /= evl;
      }
      m_div[e(0)].segment<3>(im(0) * 3) -= -ev * M_PI * m_radius(e(0) * 2 + 0) *
                                           m_radius(e(0) * 2 + 1) / m_vol(e(0));
      m_div[e(1)].segment<3>(im(1) * 3) -= ev * M_PI * m_radius(e(1) * 2 + 0) *
                                           m_radius(e(1) * 2 + 1) / m_vol(e(1));

    } else if (gidx < num_edges + num_triangles) {
      const int fidx = gidx - num_edges;

      const Vector3i& im = m_face_inv_mapping[fidx];
      const auto& f = m_faces.row(fidx);
      const int im_base_0 = (int)m_particle_to_edge[f[0]].size();
      const int im_base_1 = (int)m_particle_to_edge[f[1]].size();
      const int im_base_2 = (int)m_particle_to_edge[f[2]].size();

      Vector3s d0, d1, d2;

      mathutils::get_div_triangle(
          m_vol(f[0]), m_vol(f[1]), m_vol(f[2]),
          m_radius(f[0] * 2 + 0) + m_radius(f[0] * 2 + 1),
          m_radius(f[1] * 2 + 0) + m_radius(f[1] * 2 + 1),
          m_radius(f[2] * 2 + 0) + m_radius(f[2] * 2 + 1),
          m_x.segment<3>(f[0] * 4), m_x.segment<3>(f[1] * 4),
          m_x.segment<3>(f[2] * 4), d0, d1, d2);

      m_div[f[0]].segment<3>((im_base_0 + im(0)) * 3) += d0;
      m_div[f[1]].segment<3>((im_base_1 + im(1)) * 3) += d1;
      m_div[f[2]].segment<3>((im_base_2 + im(2)) * 3) += d2;
    }
  });
}

const std::vector<std::vector<RayTriInfo> >& TwoDScene::getIntersections()
    const {
  return m_ray_tri_gauss;
}

/*!
 * shoot rays and compute hitting points on elements. This is used for finding
 * cohesion pairs.
 */
void TwoDScene::updateIntersection() {
  const int num_edges = getNumEdges();
  const int num_faces = getNumFaces();
  const int num_soft_elasto = num_faces + num_edges;

  threadutils::for_each(0, num_soft_elasto,
    [&](int gidx) { m_ray_tri_gauss[gidx].resize(0); });
}

const std::vector<int> TwoDScene::getParticleGroup() const {
  return m_particle_group;
}

/*!
 * initialize variables on face/edge Elements
 */
void TwoDScene::initGaussSystem() {
  // update connectivity
  const int num_edges = m_edges.rows();
  const int num_triangles = m_faces.rows();
  const int num_surfels = m_surfels.size();

  // sample Gauss from edges & triangles
  const int num_system = num_edges + num_triangles + num_surfels;

  m_x_gauss.resize(4 * num_system);
  m_v_gauss.resize(4 * num_system);
  m_dv_gauss.resize(4 * num_system);
  m_fluid_v_gauss.resize(4 * num_system);
  m_m_gauss.resize(4 * num_system);
  m_fluid_m_gauss.resize(4 * num_system);
  m_vol_gauss.resize(num_system);
  m_rest_vol_gauss.resize(num_system);
  m_fluid_vol_gauss.resize(num_system);
  m_radius_gauss.resize(2 * num_system);
  m_volume_fraction_gauss.resize(num_system);
  m_rest_volume_fraction_gauss.resize(num_system);
  m_ray_tri_gauss.resize(num_system);

  m_gauss_nodes_x.resize(num_system);
  m_gauss_nodes_y.resize(num_system);
  m_gauss_nodes_z.resize(num_system);

  m_gauss_weights.resize(num_system);
  m_gauss_to_parameters.resize(num_system);

  m_Fe_gauss.resize(3 * num_system, 3);
  m_d_gauss.resize(3 * num_system, 3);
  m_d_old_gauss.resize(3 * num_system, 3);
  m_D_inv_gauss.resize(3 * num_system, 3);
  m_D_gauss.resize(3 * num_system, 3);
  m_dFe_gauss.resize(3 * num_system, 3);
  m_norm_gauss.resize(3 * num_system, 3);

  m_grad_gauss.resize(3 * num_system, 3);

  threadutils::for_each(0, num_system, [&](int i) {
    m_Fe_gauss.block<3, 3>(i * 3, 0).setIdentity();
  });

  threadutils::for_each(0, num_edges, [&](int i) {
    const auto& e = m_edges.row(i);
    m_gauss_to_parameters[i] = m_edge_to_parameters[i];
    m_x_gauss.segment<4>(i * 4) =
        (m_x.segment<4>(e(0) * 4) + m_x.segment<4>(e(1) * 4)) * 0.5;
    m_v_gauss.segment<4>(i * 4) =
        (m_v.segment<4>(e(0) * 4) + m_v.segment<4>(e(1) * 4)) * 0.5;
    m_dv_gauss.segment<4>(i * 4) =
        (m_dv.segment<4>(e(0) * 4) + m_dv.segment<4>(e(1) * 4)) * 0.5;
    m_fluid_v_gauss.segment<4>(i * 4) =
        (m_fluid_v.segment<4>(e(0) * 4) + m_fluid_v.segment<4>(e(1) * 4)) * 0.5;
    const scalar g_radius_A = getGaussRadius(i, 0);
    const scalar g_radius_B = getGaussRadius(i, 1);
    m_rest_vol_gauss(i) = m_vol_gauss(i) =
        m_edge_rest_length(i) * g_radius_A * g_radius_B * M_PI;
    m_m_gauss.segment<3>(i * 4).setConstant(m_vol_gauss(i) *
                                            getGaussDensity(i));
    m_m_gauss(i * 4 + 3) =
        m_vol_gauss(i) * getGaussDensity(i) * 0.5 * g_radius_A * g_radius_B;
    m_fluid_vol_gauss(i) = (m_fluid_vol(e(0)) + m_fluid_vol(e(1))) * 0.5;
    m_fluid_m_gauss.segment<4>(i * 4) =
        (m_fluid_m.segment<4>(e(0) * 4) + m_fluid_m.segment<4>(e(1) * 4)) * 0.5;
    m_volume_fraction_gauss(i) =
        (m_volume_fraction(e(0)) + m_volume_fraction(e(1))) * 0.5;
    m_rest_volume_fraction_gauss(i) =
        (m_rest_volume_fraction(e(0)) + m_rest_volume_fraction(e(1))) * 0.5;
    m_radius_gauss(i * 2 + 0) =
        sqrt((m_radius(e(0) * 2 + 0) * m_radius(e(0) * 2 + 0) +
              m_radius(e(1) * 2 + 0) * m_radius(e(1)) * 2 + 0) *
             0.5);
    m_radius_gauss(i * 2 + 1) =
        sqrt((m_radius(e(0) * 2 + 1) * m_radius(e(0) * 2 + 1) +
              m_radius(e(1) * 2 + 1) * m_radius(e(1)) * 2 + 1) *
             0.5);

    m_grad_gauss.block<3, 3>(i * 3, 0).setZero();
    Vector3s ev = m_x.segment<3>(e(1) * 4) - m_x.segment<3>(e(0) * 4);
    scalar l2ev = ev.squaredNorm();
    if (l2ev > 1e-20) ev /= l2ev;
    m_grad_gauss.block<3, 1>(i * 3, 0) = -ev;
    m_grad_gauss.block<3, 1>(i * 3, 1) = ev;
  });

  threadutils::for_each(0, num_triangles, [&](int i) {
    const auto& f = m_faces.row(i);
    m_gauss_to_parameters[i + num_edges] = m_face_to_parameters[i];
    const Vector3s& angle_frac = m_face_weights[i];

    m_x_gauss.segment<4>((i + num_edges) * 4) =
        m_x.segment<4>(f[0] * 4) * angle_frac[0] +
        m_x.segment<4>(f[1] * 4) * angle_frac[1] +
        m_x.segment<4>(f[2] * 4) * angle_frac[2];
    m_v_gauss.segment<4>((i + num_edges) * 4) =
        m_v.segment<4>(f[0] * 4) * angle_frac[0] +
        m_v.segment<4>(f[1] * 4) * angle_frac[1] +
        m_v.segment<4>(f[2] * 4) * angle_frac[2];
    m_dv_gauss.segment<4>((i + num_edges) * 4) =
        m_dv.segment<4>(f[0] * 4) * angle_frac[0] +
        m_dv.segment<4>(f[1] * 4) * angle_frac[1] +
        m_dv.segment<4>(f[2] * 4) * angle_frac[2];
    m_fluid_v_gauss.segment<4>((i + num_edges) * 4) =
        m_fluid_v.segment<4>(f[0] * 4) * angle_frac[0] +
        m_fluid_v.segment<4>(f[1] * 4) * angle_frac[1] +
        m_fluid_v.segment<4>(f[2] * 4) * angle_frac[2];
    const scalar g_radius_A = getGaussRadius(i + num_edges, 0);
    const scalar g_radius_B = getGaussRadius(i + num_edges, 1);
    m_rest_vol_gauss(i + num_edges) = m_vol_gauss(i + num_edges) =
        m_face_rest_area(i) * (g_radius_A + g_radius_B);
    m_m_gauss.segment<3>((i + num_edges) * 4)
        .setConstant(m_vol_gauss(i + num_edges) *
                     getGaussDensity(i + num_edges));
    m_m_gauss((i + num_edges) * 4 + 3) = 1.0;

    m_radius_gauss((i + num_edges) * 2 + 0) =
        sqrt((m_radius(f(0) * 2 + 0) * m_radius(f(0) * 2 + 0) +
              m_radius(f(1) * 2 + 0) * m_radius(f(1) * 2 + 0) +
              m_radius(f(2) * 2 + 0) * m_radius(f(2) * 2 + 0)) /
             3.0);
    m_radius_gauss((i + num_edges) * 2 + 1) =
        sqrt((m_radius(f(0) * 2 + 1) * m_radius(f(0) * 2 + 1) +
              m_radius(f(1) * 2 + 1) * m_radius(f(1) * 2 + 1) +
              m_radius(f(2) * 2 + 1) * m_radius(f(2) * 2 + 1)) /
             3.0);

    m_fluid_vol_gauss(i + num_edges) = m_fluid_vol(f[0]) * angle_frac[0] +
                                       m_fluid_vol(f[1]) * angle_frac[1] +
                                       m_fluid_vol(f[2]) * angle_frac[2];
    m_fluid_m_gauss.segment<4>((i + num_edges) * 4) =
        m_fluid_m.segment<4>(f[0] * 4) * angle_frac[0] +
        m_fluid_m.segment<4>(f[1] * 4) * angle_frac[1] +
        m_fluid_m.segment<4>(f[2] * 4) * angle_frac[2];
    m_volume_fraction_gauss(i + num_edges) =
        m_volume_fraction(f[0]) * angle_frac[0] +
        m_volume_fraction(f[1]) * angle_frac[1] +
        m_volume_fraction(f[2]) * angle_frac[2];
    m_rest_volume_fraction_gauss(i + num_edges) =
        m_rest_volume_fraction(f[0]) * angle_frac[0] +
        m_rest_volume_fraction(f[1]) * angle_frac[1] +
        m_rest_volume_fraction(f[2]) * angle_frac[2];

    Matrix3s g;
    mathutils::grad_triangle(m_rest_x.segment<3>(f[0] * 4),
                             m_rest_x.segment<3>(f[1] * 4),
                             m_rest_x.segment<3>(f[2] * 4), g);
    m_grad_gauss.block<3, 3>((i + num_edges) * 3, 0) = g;
  });

  threadutils::for_each(0, num_surfels, [&](int i) {
    int pidx = m_surfels[i];
    int gidx = i + num_edges + num_triangles;
    m_gauss_to_parameters[gidx] = 0;
    m_x_gauss.segment<4>(gidx * 4) = m_x.segment<4>(pidx * 4);
    m_v_gauss.segment<4>(gidx * 4) = m_v.segment<4>(pidx * 4);
    m_dv_gauss.segment<4>(gidx * 4) = m_dv.segment<4>(pidx * 4);
    m_fluid_v_gauss.segment<4>(gidx * 4) = m_fluid_v.segment<4>(pidx * 4);
    m_rest_vol_gauss(gidx) = m_vol_gauss(gidx) = m_vol(pidx);
    m_radius_gauss.segment<2>(gidx * 2) = m_radius.segment<2>(pidx * 2);
    m_m_gauss.segment<4>(gidx) = m_m.segment<4>(pidx);

    m_fluid_vol_gauss(gidx) = 0.0;
    m_fluid_m_gauss.segment<4>(gidx * 4).setZero();
    m_rest_volume_fraction_gauss(gidx) = m_volume_fraction_gauss(gidx) = 1.0;
    m_grad_gauss.block<3, 3>(gidx * 3, 0).setZero();
  });

  // init m_D_inv_gauss and m_D_gauss
  threadutils::for_each(0, num_edges, [&](int i) {
    const auto& e = m_edges.row(i);

    Vector3s tangent = (m_x.segment<3>(e(1) * 4) - m_x.segment<3>(e(0) * 4));
    Vector3s normal = findNormal(tangent.normalized());
    Vector3s binorm = tangent.cross(normal).normalized();

    // warning: this is a hack!!
    //        binorm = Vector3s(0, 0, 1);
    //        normal = binorm.cross(tangent).normalized();

    Matrix3s m_D;
    m_D.block<3, 1>(0, 0) = tangent;
    m_D.block<3, 1>(0, 1) = normal;
    m_D.block<3, 1>(0, 2) = binorm;

    m_d_gauss.block<3, 3>(i * 3, 0) = m_D;

    m_norm_gauss.block<3, 1>(0, 0) = tangent.normalized();
    m_norm_gauss.block<3, 1>(0, 1) = normal.normalized();
    m_norm_gauss.block<3, 1>(0, 2) = binorm.normalized();

    Matrix3s Dstar = Matrix3s::Identity();
    Dstar(0, 0) = tangent.norm();
    Matrix3s invDstar = Dstar.inverse();

    m_Fe_gauss.block<3, 3>(i * 3, 0) = m_D * invDstar;

    m_D_inv_gauss.block<3, 3>(i * 3, 0) = invDstar;

    m_D_gauss.block<3, 3>(i * 3, 0) = Dstar;
  });

  threadutils::for_each(num_edges, num_edges + num_triangles, [&](int i) {
    const auto& f = m_faces.row(i - num_edges);
    Vector3s t0 = (m_x.segment<3>(f[1] * 4) - m_x.segment<3>(f[0] * 4));
    Vector3s t1 = (m_x.segment<3>(f[2] * 4) - m_x.segment<3>(f[0] * 4));
    Vector3s norm = t1.cross(t0).normalized();

    Matrix3s m_D;
    m_D.block<3, 1>(0, 0) = t0;
    m_D.block<3, 1>(0, 1) = t1;
    m_D.block<3, 1>(0, 2) = norm;

    m_d_gauss.block<3, 3>(i * 3, 0) = m_D;

    Matrix3s Q, R;

    mathutils::QRDecompose<scalar, 3>(m_D, Q, R);

    m_norm_gauss.block<3, 1>(i * 3, 0) = t0.normalized();
    m_norm_gauss.block<3, 1>(i * 3, 1) = t0.cross(norm).normalized();
    m_norm_gauss.block<3, 1>(i * 3, 2) = norm;

    // compute rotations - from norm to Z
    Eigen::Quaternion<scalar> rot0 =
        Eigen::Quaternion<scalar>::FromTwoVectors(norm, Vector3s::UnitZ());

    Vector3s u = rot0 * t0;
    Vector3s v = rot0 * t1;

    // compute rotations - from t0 to X
    Eigen::Quaternion<scalar> rot1 =
        Eigen::Quaternion<scalar>::FromTwoVectors(u, Vector3s::UnitX());

    Vector3s ru = rot1 * u;
    Vector3s rv = rot1 * v;

    Matrix3s Dstar = Matrix3s::Identity();
    Dstar(0, 0) = ru(0);
    Dstar(0, 1) = rv(0);
    Dstar(1, 1) = rv(1);

    Matrix3s invDstar = Dstar.inverse();

    m_Fe_gauss.block<3, 3>(i * 3, 0) = m_D * invDstar;

    m_D_inv_gauss.block<3, 3>(i * 3, 0) = invDstar;

    m_D_gauss.block<3, 3>(i * 3, 0) = Dstar;
  });

  threadutils::for_each(num_edges + num_triangles,
                        num_edges + num_triangles + num_surfels, [&](int i) {
                          const int s = i - num_edges - num_triangles;
                          const Vector3s& norm = m_surfel_norms[s];
                          Eigen::Quaternion<scalar> rot0 =
                              Eigen::Quaternion<scalar>::FromTwoVectors(
                                  Vector3s::UnitZ(), norm);

                          Matrix3s m_D;
                          m_D.block<3, 1>(0, 0) = rot0 * Vector3s::UnitX();
                          m_D.block<3, 1>(0, 1) = rot0 * Vector3s::UnitY();
                          m_D.block<3, 1>(0, 2) = norm;

                          m_norm_gauss.block<3, 3>(i * 3, 0) = m_D;

                          m_d_gauss.block<3, 3>(i * 3, 0) = m_D;
                          m_Fe_gauss.block<3, 3>(i * 3, 0) = m_D;
                          m_D_inv_gauss.block<3, 3>(i * 3, 0).setIdentity();
                          m_D_gauss.block<3, 3>(i * 3, 0).setIdentity();
                        });
}

/*!
 * compute derivative of energy E over deformation gradient Fe, this is crucial
 * for computing collision force
 */
void TwoDScene::computedEdFe() {
  const int num_gauss = getNumGausses();
  const int num_edges = getNumEdges();

  m_dFe_gauss.setZero();

  // compute forces on yarns
  threadutils::for_each(0, num_edges, [&](int i) {
    Matrix3s FeD = m_d_gauss.block<3, 3>(i * 3, 0);
    Matrix3s Q, R;

    mathutils::QRDecompose<scalar, 3>(FeD, Q, R);

    double dhdr22, dhdr23, dhdr33;
    double mu, la;
    mu = getMu(i) * getCollisionMultiplier(i);
    la = getLa(i) * getCollisionMultiplier(i);

    mathutils::dhdr_yarn(mu, la, R(1, 1), R(1, 2), R(2, 2), &dhdr22, &dhdr23,
                         &dhdr33);

    // dphidr = [f' g'rt]
    //         [0  P_hat]

    Matrix3s dphidr;
    dphidr.setZero();

    if (dhdr22 != 0.0 || dhdr33 != 0.0) {
      dphidr(0, 1) = mu * R(0, 1);
      dphidr(0, 2) = mu * R(0, 2);
    }

    dphidr(1, 1) = dhdr22;
    dphidr(2, 2) = dhdr33;
    dphidr(1, 2) = dhdr23;

    const Matrix3s& dphidrRt = dphidr * R.transpose();

    const Matrix3s& tauK = dphidrRt.triangularView<Eigen::Upper>();
    const Matrix3s& tauKt = tauK.transpose();
    const Matrix3s& DK = dphidrRt.diagonal().asDiagonal();
    const Matrix3s& dphidd = Q * (tauK + tauKt - DK) * R.inverse().transpose();

    assert(!std::isnan(dphidd.sum()));
    assert(!std::isinf(dphidd.sum()));

    m_dFe_gauss.block<3, 3>(i * 3, 0) =
        dphidd * m_D_gauss.block<3, 3>(i * 3, 0).transpose();
  });

  // compute forces on clothes and surfels (removed)
  
}

/*!
 * update bounding box
 */
void TwoDScene::updateParticleBoundingBox() {
  Vector4s bbmin = Vector4s::Constant(1e+20);
  Vector4s bbmax = Vector4s::Constant(-1e+20);

  bbmin = threadutils::reduction((Vector4s*)m_x.data(), getNumParticles(),
                                 bbmin, [](Vector4s x, Vector4s y) -> Vector4s {
                                   return Vector4s(std::min(x(0), y(0)),
                                                   std::min(x(1), y(1)),
                                                   std::min(x(2), y(2)), 0.0);
                                 });

  bbmax = threadutils::reduction((Vector4s*)m_x.data(), getNumParticles(),
                                 bbmax, [](Vector4s x, Vector4s y) -> Vector4s {
                                   return Vector4s(std::max(x(0), y(0)),
                                                   std::max(x(1), y(1)),
                                                   std::max(x(2), y(2)), 0.0);
                                 });

  const scalar dx = m_bucket_size * 2.0;

  m_bbx_min = Vector3s(floor(bbmin(0) / dx) * dx, floor(bbmin(1) / dx) * dx,
                       floor(bbmin(2) / dx) * dx);
  m_bbx_max = Vector3s(ceil(bbmax(0) / dx) * dx, ceil(bbmax(1) / dx) * dx,
                       ceil(bbmax(2) / dx) * dx);
}

void TwoDScene::setBucketInfo(const scalar& bucket_size, int num_nodes,
                              int kernel_order) {
  m_bucket_size = bucket_size;
  m_num_nodes = num_nodes;
  m_num_bucket_colors = 2;
}

int TwoDScene::getNumColors() const { return m_num_colors; }

const std::vector<VectorXi>& TwoDScene::getNodeColorP() const {
  return m_node_color_p;
}

const std::vector<std::pair<int, int> >& TwoDScene::getNodeParticlePairsX(
    int bucket_idx, int pidx) const {
  return m_node_particles_x[bucket_idx][pidx];
}

const std::vector<std::pair<int, int> >& TwoDScene::getNodeParticlePairsY(
    int bucket_idx, int pidx) const {
  return m_node_particles_y[bucket_idx][pidx];
}

const std::vector<std::pair<int, int> >& TwoDScene::getNodeParticlePairsZ(
    int bucket_idx, int pidx) const {
  return m_node_particles_z[bucket_idx][pidx];
}

int TwoDScene::getNumBucketColors() const { return m_num_bucket_colors; }

int TwoDScene::getKernelOrder() const { return m_kernel_order; }

scalar TwoDScene::getBucketLength() const {
  return getCellSize() * (scalar)m_num_nodes;
}

const Vector3s& TwoDScene::getBucketMinCorner() const {
  return m_bucket_mincorner;
}

/*!
 * Put particles into buckets for parallel searching and computing
 */
void TwoDScene::rebucketizeParticles() {
  scalar dx = getCellSize();

  const scalar extra_border = 3.0;

  Vector3s content_size =
      m_bbx_max - m_bbx_min +
      Vector3s::Constant(m_bucket_size * extra_border * 2.0);

  Vector3i grid_num_cells =
      Vector3i(std::max(1, (int)ceil(content_size(0) / dx)),
               std::max(1, (int)ceil(content_size(1) / dx)),
               std::max(1, (int)ceil(content_size(2) / dx)));

  Vector3s grid_size =
      Vector3s((scalar)grid_num_cells[0] * dx, (scalar)grid_num_cells[1] * dx,
               (scalar)grid_num_cells[2] * dx);

  m_grid_mincorner =
      m_bbx_min - Vector3s::Constant(m_bucket_size * extra_border);
  m_bucket_mincorner = m_grid_mincorner;

  Vector3i num_buckets =
      Vector3i(std::max(1, (int)ceil(grid_size(0) / m_bucket_size)),
               std::max(1, (int)ceil(grid_size(1) / m_bucket_size)),
               std::max(1, (int)ceil(grid_size(2) / m_bucket_size)));

  m_particle_buckets.resize(num_buckets(0), num_buckets(1), num_buckets(2));
  m_gauss_buckets.resize(num_buckets(0), num_buckets(1), num_buckets(2));
  m_particle_cells.resize(grid_num_cells(0), grid_num_cells(1),
                          grid_num_cells(2));

  m_particle_buckets.sort(getNumParticles(), [&](int pidx, int& i, int& j,
                                                 int& k) {
    i = (int)floor((m_x(pidx * 4 + 0) - m_bucket_mincorner(0)) / m_bucket_size);
    j = (int)floor((m_x(pidx * 4 + 1) - m_bucket_mincorner(1)) / m_bucket_size);
    k = (int)floor((m_x(pidx * 4 + 2) - m_bucket_mincorner(2)) / m_bucket_size);
  });

  m_gauss_buckets.sort(getNumGausses(), [&](int pidx, int& i, int& j, int& k) {
    i = (int)floor((m_x_gauss(pidx * 4 + 0) - m_bucket_mincorner(0)) /
                   m_bucket_size);
    j = (int)floor((m_x_gauss(pidx * 4 + 1) - m_bucket_mincorner(1)) /
                   m_bucket_size);
    k = (int)floor((m_x_gauss(pidx * 4 + 2) - m_bucket_mincorner(2)) /
                   m_bucket_size);
  });

  const int total_buckets = m_particle_buckets.size();

  m_bucket_activated.assign(total_buckets, 0U);
}

/*!
 * remove empty particles.
 */
void TwoDScene::removeEmptyParticles() {
  const int num_parts = getNumParticles();
  const int num_elasto = getNumElastoParticles();

  int new_num_parts = num_parts;
  for (int i = num_elasto; i < new_num_parts;) {
    if (m_fluid_vol(i) < 1e-20) {
      swapParticles(i, --new_num_parts);
    } else {
      ++i;
    }
  }

  if (new_num_parts < num_parts) {
    conservativeResizeParticles(new_num_parts);

    m_fluids.resize(new_num_parts - num_elasto);
    for (int i = num_elasto; i < new_num_parts; ++i) {
      m_fluids[i - num_elasto] = i;
    }

    m_particle_buckets.sort(
        new_num_parts, [&](int pidx, int& i, int& j, int& k) {
          i = (int)floor((m_x(pidx * 4 + 0) - m_bucket_mincorner(0)) /
                         m_bucket_size);
          j = (int)floor((m_x(pidx * 4 + 1) - m_bucket_mincorner(1)) /
                         m_bucket_size);
          k = (int)floor((m_x(pidx * 4 + 2) - m_bucket_mincorner(2)) /
                         m_bucket_size);
        });
  }
}

void TwoDScene::updateStrandParamViscosity(const scalar& dt) {
  const int num_params = (int)m_strandParameters.size();
  threadutils::for_each(0, num_params, [&](int i) {
    m_strandParameters[i]->computeViscousForceCoefficients(dt);
  });
}

/*!
 * when particles are outside a Terminator, we delete them
 */
void TwoDScene::terminateParticles() {

  const int num_parts = getNumParticles();
  const int num_elasto = getNumElastoParticles();
  threadutils::for_each(num_elasto, num_parts, [&](int pidx) {
    const Vector3s& pos = m_x.segment<3>(pidx * 4);
    Vector3s vel;
    const scalar phi = computePhiVel(pos, vel);
    if (phi < 0.0) m_fluid_vol(pidx) = 0.0;
  });

  removeEmptyParticles();
}

/*!
 * project particles to avoid penetrating rigid bodies
 */
void TwoDScene::solidProjection(const scalar& dt) {
  const int num_parts = getNumParticles();
  const int num_elasto = getNumElastoParticles();

  const scalar iD = getInverseDCoeff();
  // solid projection 
  threadutils::for_each(num_elasto, num_parts, [&](int pidx) {
    if (m_particle_to_surfel[pidx] >= 0) return;

    const auto& node_indices_sphi = m_particle_nodes_solid_phi[pidx];
    const auto& particle_weights = m_particle_weights[pidx];

    scalar phi_ori = 0.0;
    Vector3s grad_phi = Vector3s::Zero();
    const Vector3s& pos = m_x.segment<3>(pidx * 4);

    for (int nidx = 0; nidx < node_indices_sphi.rows(); ++nidx) {
      const int bucket_idx = node_indices_sphi(nidx, 0);
      const int node_idx = node_indices_sphi(nidx, 1);

      scalar phi;
      if (m_bucket_activated[bucket_idx]) {
        phi = m_node_solid_phi[bucket_idx](node_idx);
      } else {
        phi = 3.0 * getCellSize();
      }

      const scalar w = particle_weights(nidx, 3);
      const Vector3s& np = m_node_pos[bucket_idx].segment<3>(node_idx * 3);

      phi_ori += phi * w;
      grad_phi += phi * iD * w * (np - pos);
    }

    if (grad_phi.norm() > 1e-20) grad_phi.normalize();

    const Vector3s dpos = m_fluid_v.segment<3>(pidx * 4) * dt;

    scalar phi_now = phi_ori + grad_phi.dot(dpos);

    //            if(pidx < num_elasto) phi_now *= 0.1;

    if (phi_now < 0.0) {
      m_x.segment<3>(pidx * 4) -= phi_now * grad_phi;
    }
  });
}

/*!
 * compute particle weight on solid node
 */
void TwoDScene::updateSolidWeights() {
  const int num_buckets = m_particle_buckets.size();
  m_node_solid_weight_x.resize(num_buckets);
  m_node_solid_weight_y.resize(num_buckets);
  m_node_solid_weight_z.resize(num_buckets);

  const scalar dx = getCellSize();

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    const VectorXi& bucket_node_idx_solid_phi_x =
        m_node_index_solid_phi_x[bucket_idx];
    const VectorXi& bucket_node_idx_solid_phi_y =
        m_node_index_solid_phi_y[bucket_idx];
    const VectorXi& bucket_node_idx_solid_phi_z =
        m_node_index_solid_phi_z[bucket_idx];

    const int num_solid_phi = getNumNodes(bucket_idx);

    VectorXs& bucket_weight_x = m_node_solid_weight_x[bucket_idx];
    VectorXs& bucket_weight_y = m_node_solid_weight_y[bucket_idx];
    VectorXs& bucket_weight_z = m_node_solid_weight_z[bucket_idx];

    if (bucket_weight_x.size() != num_solid_phi)
      bucket_weight_x.resize(num_solid_phi);
    if (bucket_weight_y.size() != num_solid_phi)
      bucket_weight_y.resize(num_solid_phi);
    if (bucket_weight_z.size() != num_solid_phi)
      bucket_weight_z.resize(num_solid_phi);

    for (int i = 0; i < num_solid_phi; ++i) {
      const Vector8i& indices = bucket_node_idx_solid_phi_x.segment<8>(i * 8);
      scalar phi0 = 0.5 * dx;
      scalar phi1 = 0.5 * dx;
      scalar phi2 = 0.5 * dx;
      scalar phi3 = 0.5 * dx;

      if (indices[0] >= 0 && m_bucket_activated[indices[0]])
        phi0 = m_node_solid_phi[indices[0]][indices[1]];
      if (indices[2] >= 0 && m_bucket_activated[indices[2]])
        phi1 = m_node_solid_phi[indices[2]][indices[3]];
      if (indices[4] >= 0 && m_bucket_activated[indices[4]])
        phi2 = m_node_solid_phi[indices[4]][indices[5]];
      if (indices[6] >= 0 && m_bucket_activated[indices[6]])
        phi3 = m_node_solid_phi[indices[6]][indices[7]];

      bucket_weight_x(i) = mathutils::clamp(
          1.0 - mathutils::fraction_inside(phi0, phi1, phi2, phi3), 0.0, 1.0);
    }

    for (int i = 0; i < num_solid_phi; ++i) {
      const Vector8i& indices = bucket_node_idx_solid_phi_y.segment<8>(i * 8);
      scalar phi0 = 0.5 * dx;
      scalar phi1 = 0.5 * dx;
      scalar phi2 = 0.5 * dx;
      scalar phi3 = 0.5 * dx;

      if (indices[0] >= 0 && m_bucket_activated[indices[0]])
        phi0 = m_node_solid_phi[indices[0]][indices[1]];
      if (indices[2] >= 0 && m_bucket_activated[indices[2]])
        phi1 = m_node_solid_phi[indices[2]][indices[3]];
      if (indices[4] >= 0 && m_bucket_activated[indices[4]])
        phi2 = m_node_solid_phi[indices[4]][indices[5]];
      if (indices[6] >= 0 && m_bucket_activated[indices[6]])
        phi3 = m_node_solid_phi[indices[6]][indices[7]];

      bucket_weight_y(i) = mathutils::clamp(
          1.0 - mathutils::fraction_inside(phi0, phi1, phi2, phi3), 0.0, 1.0);
    }

    for (int i = 0; i < num_solid_phi; ++i) {
      const Vector8i& indices = bucket_node_idx_solid_phi_z.segment<8>(i * 8);
      scalar phi0 = 0.5 * dx;
      scalar phi1 = 0.5 * dx;
      scalar phi2 = 0.5 * dx;
      scalar phi3 = 0.5 * dx;

      if (indices[0] >= 0 && m_bucket_activated[indices[0]])
        phi0 = m_node_solid_phi[indices[0]][indices[1]];
      if (indices[2] >= 0 && m_bucket_activated[indices[2]])
        phi1 = m_node_solid_phi[indices[2]][indices[3]];
      if (indices[4] >= 0 && m_bucket_activated[indices[4]])
        phi2 = m_node_solid_phi[indices[4]][indices[5]];
      if (indices[6] >= 0 && m_bucket_activated[indices[6]])
        phi3 = m_node_solid_phi[indices[6]][indices[7]];

      bucket_weight_z(i) = mathutils::clamp(
          1.0 - mathutils::fraction_inside(phi0, phi1, phi2, phi3), 0.0, 1.0);
    }
  });
}

Vector3s TwoDScene::nodePosFromBucket(int bucket_idx, int raw_node_idx,
                                      const Vector3s& offset) const {
  Vector3i handle = m_particle_buckets.bucket_handle(bucket_idx);
  Vector3s bucket_left_corner =
      m_bucket_mincorner + Vector3s(handle(0) * m_bucket_size,
                                    handle(1) * m_bucket_size,
                                    handle(2) * m_bucket_size);
  int iz = raw_node_idx / (m_num_nodes * m_num_nodes);
  int ixy = raw_node_idx - iz * m_num_nodes * m_num_nodes;
  int iy = ixy / m_num_nodes;
  int ix = ixy - iy * m_num_nodes;
  Vector3s node_pos =
      bucket_left_corner + (Vector3s(ix, iy, iz) + offset) * getCellSize();

  return node_pos;
}

/*!
 * update particle weights on X-, Y-, Z-, SOLID- and PRESSURE- nodes.
 */
void TwoDScene::updateParticleWeights(scalar dt, int start, int end) {
  const scalar h = getCellSize();

  threadutils::for_each(start, end, [&](int pidx) {
    if (m_inside[pidx] == 0) return;

    const Matrix27x2i& indices_x = m_particle_nodes_x[pidx];
    const Matrix27x2i& indices_y = m_particle_nodes_y[pidx];
    const Matrix27x2i& indices_z = m_particle_nodes_z[pidx];
    const Matrix27x2i& indices_sphi = m_particle_nodes_solid_phi[pidx];
    const Matrix27x2i& indices_p = m_particle_nodes_p[pidx];

    Vector27s& weights_p = m_particle_weights_p[pidx];
    Matrix27x4s& weights = m_particle_weights[pidx];

    const Vector3s& pos = m_x.segment<3>(pidx * 4);

    for (int nidx = 0; nidx < indices_p.rows(); ++nidx) {
      const int node_bucket_idx = indices_p(nidx, 0);
      const int node_idx = indices_p(nidx, 1);
      Vector3s dx;

      const Vector3s& np = getNodePosP(node_bucket_idx, node_idx);
      dx = (pos - np) / h;

      const scalar w = mathutils::N_kernel<2>(dx);

      weights_p(nidx) = w;
    }

    for (int nidx = 0; nidx < indices_sphi.rows(); ++nidx) {
      const int node_bucket_idx = indices_sphi(nidx, 0);
      const int node_idx = indices_sphi(nidx, 1);

      Vector3s dx;

      const Vector3s& np = getNodePosSolidPhi(node_bucket_idx, node_idx);
      dx = (pos - np) / h;

      Vector3s g;
      const scalar w = mathutils::N_kernel<2>(dx);

      weights(nidx, 3) = w;
    }

    for (int nidx = 0; nidx < indices_x.rows(); ++nidx) {
      const int node_bucket_idx = indices_x(nidx, 0);
      const int node_idx = indices_x(nidx, 1);

      Vector3s dx;

      const Vector3s& np = getNodePosX(node_bucket_idx, node_idx);
      dx = (pos - np) / h;

      Vector3s g;
      const scalar w = mathutils::N_kernel<2>(dx);

      weights(nidx, 0) = w;
    }

    for (int nidx = 0; nidx < indices_y.rows(); ++nidx) {
      const int node_bucket_idx = indices_y(nidx, 0);
      const int node_idx = indices_y(nidx, 1);

      Vector3s dx;

      const Vector3s& np = getNodePosY(node_bucket_idx, node_idx);
      dx = (pos - np) / h;

      Vector3s g;
      const scalar w = mathutils::N_kernel<2>(dx);

      weights(nidx, 1) = w;
    }

    for (int nidx = 0; nidx < indices_z.rows(); ++nidx) {
      const int node_bucket_idx = indices_z(nidx, 0);
      const int node_idx = indices_z(nidx, 1);

      Vector3s dx;

      const Vector3s& np = getNodePosZ(node_bucket_idx, node_idx);
      dx = (pos - np) / h;

      Vector3s g;
      const scalar w = mathutils::N_kernel<2>(dx);

      weights(nidx, 2) = w;
    }
  });
}

/*!
 * update element weight on nodes
 */
void TwoDScene::updateGaussWeights(scalar dt) {
  const int num_gauss = getNumGausses();
  const int num_edges = getNumEdges();
  const int num_faces = getNumFaces();
  const int num_soft_gauss = num_edges + num_faces;
  const scalar h = getCellSize();

  threadutils::for_each(0, num_gauss, [&](int pidx) {
    if (pidx >= num_soft_gauss &&
        m_inside[m_surfels[pidx - num_soft_gauss]] == 0)
      return;

    const Matrix27x2i& indices_x = m_gauss_nodes_x[pidx];
    const Matrix27x2i& indices_y = m_gauss_nodes_y[pidx];
    const Matrix27x2i& indices_z = m_gauss_nodes_z[pidx];

    Matrix27x3s& weights = m_gauss_weights[pidx];

    const Vector3s& pos = m_x_gauss.segment<3>(pidx * 4);

    for (int nidx = 0; nidx < indices_x.rows(); ++nidx) {
      const int node_bucket_idx = indices_x(nidx, 0);
      const int node_idx = indices_x(nidx, 1);

      Vector3s dx;

      const Vector3s& np = getNodePosX(node_bucket_idx, node_idx);
      dx = (pos - np) / h;

      Vector3s g;
      const scalar w = mathutils::N_kernel<2>(dx);

      weights(nidx, 0) = w;
    }

    for (int nidx = 0; nidx < indices_y.rows(); ++nidx) {
      const int node_bucket_idx = indices_y(nidx, 0);
      const int node_idx = indices_y(nidx, 1);

      Vector3s dx;

      const Vector3s& np = getNodePosY(node_bucket_idx, node_idx);
      dx = (pos - np) / h;

      Vector3s g;
      const scalar w = mathutils::N_kernel<2>(dx);

      weights(nidx, 1) = w;
    }

    for (int nidx = 0; nidx < indices_z.rows(); ++nidx) {
      const int node_bucket_idx = indices_z(nidx, 0);
      const int node_idx = indices_z(nidx, 1);

      Vector3s dx;

      const Vector3s& np = getNodePosZ(node_bucket_idx, node_idx);
      dx = (pos - np) / h;

      Vector3s g;
      const scalar w = mathutils::N_kernel<2>(dx);

      weights(nidx, 2) = w;
    }
  });
}

void TwoDScene::computeWeights(scalar dt) {
  updateParticleWeights(dt, 0, getNumParticles());

  updateGaussWeights(dt);

  buildNodeParticlePairs();
}

/*!
 * from nodes, find neighbor particles, and record them
 */
void TwoDScene::buildNodeParticlePairs() {
  const int num_buckets = (int)m_particle_buckets.size();

  if ((int)m_node_particles_x.size() != num_buckets)
    m_node_particles_x.resize(num_buckets);
  if ((int)m_node_particles_y.size() != num_buckets)
    m_node_particles_y.resize(num_buckets);
  if ((int)m_node_particles_z.size() != num_buckets)
    m_node_particles_z.resize(num_buckets);
  if ((int)m_node_particles_p.size() != num_buckets)
    m_node_particles_p.resize(num_buckets);

  // re-allocate space
  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    int num_nodes = getNumNodes(bucket_idx);

    auto& bucket_node_particles_x = m_node_particles_x[bucket_idx];
    if ((int)bucket_node_particles_x.size() != num_nodes)
      bucket_node_particles_x.resize(num_nodes);

    for (int i = 0; i < num_nodes; ++i) {
      int osize = bucket_node_particles_x[i].size();
      bucket_node_particles_x[i].resize(0);
      bucket_node_particles_x[i].reserve(osize);
    }

    auto& bucket_node_particles_y = m_node_particles_y[bucket_idx];
    if ((int)bucket_node_particles_y.size() != num_nodes)
      bucket_node_particles_y.resize(num_nodes);

    for (int i = 0; i < num_nodes; ++i) {
      int osize = bucket_node_particles_y[i].size();
      bucket_node_particles_y[i].resize(0);
      bucket_node_particles_y[i].reserve(osize);
    }

    auto& bucket_node_particles_z = m_node_particles_z[bucket_idx];
    if ((int)bucket_node_particles_z.size() != num_nodes)
      bucket_node_particles_z.resize(num_nodes);

    for (int i = 0; i < num_nodes; ++i) {
      int osize = bucket_node_particles_z[i].size();
      bucket_node_particles_z[i].resize(0);
      bucket_node_particles_z[i].reserve(osize);
    }

    auto& bucket_node_particles_p = m_node_particles_p[bucket_idx];
    if ((int)bucket_node_particles_p.size() != num_nodes)
      bucket_node_particles_p.resize(num_nodes);

    for (int i = 0; i < num_nodes; ++i) {
      int osize = bucket_node_particles_p[i].size();
      bucket_node_particles_p[i].resize(0);
      bucket_node_particles_p[i].reserve(osize);
    }
  });

  m_particle_buckets.for_each_bucket_particles_colored(
      [&](int pidx, int bucket_idx) {
        if (!m_bucket_activated[bucket_idx]) return;

        auto& indices_x = getParticleNodesX(pidx);
        auto& indices_y = getParticleNodesY(pidx);
        auto& indices_z = getParticleNodesZ(pidx);
        auto& indices_p = m_particle_nodes_p[pidx];

        auto& weights = m_particle_weights[pidx];
        auto& weights_p = m_particle_weights_p[pidx];

        for (int i = 0; i < indices_x.rows(); ++i) {
          if (m_bucket_activated[indices_x(i, 0)] && weights(i, 0) > 0.0) {
            m_node_particles_x[indices_x(i, 0)][indices_x(i, 1)].emplace_back(
                std::pair<int, int>(pidx, i));
          }
        }

        for (int i = 0; i < indices_y.rows(); ++i) {
          if (m_bucket_activated[indices_y(i, 0)] && weights(i, 1) > 0.0) {
            m_node_particles_y[indices_y(i, 0)][indices_y(i, 1)].emplace_back(
                std::pair<int, int>(pidx, i));
          }
        }

        for (int i = 0; i < indices_z.rows(); ++i) {
          if (m_bucket_activated[indices_z(i, 0)] && weights(i, 2) > 0.0) {
            m_node_particles_z[indices_z(i, 0)][indices_z(i, 1)].emplace_back(
                std::pair<int, int>(pidx, i));
          }
        }

        for (int i = 0; i < indices_p.rows(); ++i) {
          if (m_bucket_activated[indices_p(i, 0)] && weights_p(i) > 0.0) {
            m_node_particles_p[indices_p(i, 0)][indices_p(i, 1)].emplace_back(
                std::pair<int, int>(pidx, i));
          }
        }
      },
      3);
}

const std::vector<VectorXs>& TwoDScene::getNodeSurfTensionP() const {
  return m_node_surf_tension;
}

scalar TwoDScene::interpolateValue(const Vector3s& pos,
                                   const std::vector<VectorXs>& phi,
                                   const Vector3s& phi_ori,
                                   const scalar& default_val) {
  const scalar dx = getCellSize();
  Vector3s grid_pos = pos - phi_ori;
  Vector3s base_pos = grid_pos / dx;
  Vector3i base_idx = Vector3i((int)floor(base_pos(0)), (int)floor(base_pos(1)),
                               (int)floor(base_pos(2)));

  scalar buf[8];
  for (int t = 0; t < 2; ++t)
    for (int s = 0; s < 2; ++s)
      for (int r = 0; r < 2; ++r) {
        int local_idx = t * 4 + s * 2 + r;
        Vector3i query_idx = base_idx + Vector3i(r, s, t);
        Vector3i bucket_handle =
            Vector3i(query_idx(0) / m_num_nodes, query_idx(1) / m_num_nodes,
                     query_idx(2) / m_num_nodes);

        if (bucket_handle(0) < 0 || bucket_handle(0) >= m_particle_buckets.ni ||
            bucket_handle(1) < 0 || bucket_handle(1) >= m_particle_buckets.nj ||
            bucket_handle(2) < 0 || bucket_handle(2) >= m_particle_buckets.nk) {
          buf[local_idx] = default_val;
          continue;
        }

        const int bucket_idx = m_particle_buckets.bucket_index(bucket_handle);
        if (!m_bucket_activated[bucket_idx]) {
          buf[local_idx] = default_val;
          continue;
        }

        Vector3i node_handle =
            Vector3i(query_idx(0) - bucket_handle(0) * m_num_nodes,
                     query_idx(1) - bucket_handle(1) * m_num_nodes,
                     query_idx(2) - bucket_handle(2) * m_num_nodes);

        const int node_idx = node_handle(2) * m_num_nodes * m_num_nodes +
                             node_handle(1) * m_num_nodes + node_handle(0);

        buf[local_idx] = phi[bucket_idx][node_idx];
      }

  Vector3s frac = Vector3s(base_pos(0) - (scalar)base_idx(0),
                           base_pos(1) - (scalar)base_idx(1),
                           base_pos(2) - (scalar)base_idx(2));

  return mathutils::trilerp(buf[0], buf[1], buf[2], buf[3], buf[4], buf[5],
                            buf[6], buf[7], frac[0], frac[1], frac[2]);
}

const std::vector<VectorXuc>& TwoDScene::getNodeLiquidValidX() const {
  return m_node_liquid_valid_x;
}

std::vector<VectorXuc>& TwoDScene::getNodeLiquidValidX() {
  return m_node_liquid_valid_x;
}

const std::vector<VectorXuc>& TwoDScene::getNodeLiquidValidY() const {
  return m_node_liquid_valid_y;
}

std::vector<VectorXuc>& TwoDScene::getNodeLiquidValidY() {
  return m_node_liquid_valid_y;
}

const std::vector<VectorXuc>& TwoDScene::getNodeLiquidValidZ() const {
  return m_node_liquid_valid_z;
}

std::vector<VectorXuc>& TwoDScene::getNodeLiquidValidZ() {
  return m_node_liquid_valid_z;
}

void TwoDScene::preAllocateNodes() {
  // mark buckets as active
  const int num_buckets = m_particle_buckets.size();

  m_node_pos.resize(num_buckets);

  m_node_pressure_neighbors.resize(num_buckets);
  m_node_pp_neighbors.resize(num_buckets);
  m_node_index_pressure_x.resize(num_buckets);
  m_node_index_pressure_y.resize(num_buckets);
  m_node_index_pressure_z.resize(num_buckets);

  m_node_index_solid_phi_x.resize(num_buckets);
  m_node_index_solid_phi_y.resize(num_buckets);
  m_node_index_solid_phi_z.resize(num_buckets);
}

/*!
 * for all particles, find the neighbor nodes to construct node structure
 */
template <typename Callable>
void TwoDScene::findNodes(const Sorter& buckets, const VectorXs& x,
                          std::vector<Matrix27x2i>& particle_nodes,
                          const Vector3s& offset, Callable func) {
  const scalar dx = getCellSize();

  // make connection for particles
  buckets.for_each_bucket_colored([&](int bucket_idx) {
    Vector3i bucket_handle = buckets.bucket_handle(bucket_idx);

    Vector3s cell_local_corner =
        Vector3s((scalar)bucket_handle(0) * m_bucket_size,
                 (scalar)bucket_handle(1) * m_bucket_size,
                 (scalar)bucket_handle(2) * m_bucket_size) +
        m_grid_mincorner + offset * dx;

    buckets.get_bucket(bucket_idx, [&](int pidx) {
      Matrix27x2i& indices = particle_nodes[pidx];

      Vector3s local_pos = (x.segment<3>(pidx * 4) - cell_local_corner) / dx;
      Vector3i ilocal_pos =
          Vector3i((int)floor(local_pos(0)), (int)floor(local_pos(1)),
                   (int)floor(local_pos(2)));

      Vector3s local_frac = mathutils::frac<scalar, 3, 1>(local_pos);

      const int klow = local_frac(2) > 0.5 ? 0 : -1;
      const int jlow = local_frac(1) > 0.5 ? 0 : -1;
      const int ilow = local_frac(0) > 0.5 ? 0 : -1;
      const int khigh = klow + 2;
      const int jhigh = jlow + 2;
      const int ihigh = ilow + 2;

      bool hasNewNode = func(pidx);

      for (int k = klow; k <= khigh; ++k)
        for (int j = jlow; j <= jhigh; ++j)
          for (int i = ilow; i <= ihigh; ++i) {
            Vector3i cell_local_idx = ilocal_pos + Vector3i(i, j, k);
            Vector3i node_bucket_handle = bucket_handle;

            for (int r = 0; r < 3; ++r) {
              while (cell_local_idx(r) < 0) {
                node_bucket_handle(r)--;
                cell_local_idx(r) += m_num_nodes;
              }
              while (cell_local_idx(r) >= m_num_nodes) {
                node_bucket_handle(r)++;
                cell_local_idx(r) -= m_num_nodes;
              }

              assert(cell_local_idx(r) >= 0 && cell_local_idx(r) < m_num_nodes);
              assert(node_bucket_handle(r) >= 0 &&
                     node_bucket_handle(r) < m_particle_buckets.dim_size(r));
            }

            int node_bucket_idx = buckets.bucket_index(node_bucket_handle);
            if (hasNewNode) {
              m_bucket_activated[node_bucket_idx] = 1U;
            }

            int cell_idx = cell_local_idx(2) * m_num_nodes * m_num_nodes +
                           cell_local_idx(1) * m_num_nodes + cell_local_idx(0);

            int nidx = (k - klow) * (3 * 3) + (j - jlow) * 3 + (i - ilow);
            indices(nidx, 0) = node_bucket_idx;
            indices(nidx, 1) = cell_idx;
          }
    });
  });
}

Vector3i TwoDScene::getNodeHandle(int node_idx) const {
  int iz = node_idx / (m_num_nodes * m_num_nodes);
  node_idx -= iz * m_num_nodes * m_num_nodes;
  int iy = node_idx / m_num_nodes;
  int ix = node_idx - iy * m_num_nodes;

  return Vector3i(ix, iy, iz);
}

int TwoDScene::getNodeIndex(const Vector3i& handle) const {
  return handle(2) * m_num_nodes * m_num_nodes + handle(1) * m_num_nodes +
         m_num_nodes;
}

Vector3s TwoDScene::getNodePosSolidPhi(int bucket_idx, int node_idx) const {
  if (m_bucket_activated[bucket_idx])
    return m_node_pos[bucket_idx].segment<3>(node_idx * 3);
  else
    return nodePosFromBucket(bucket_idx, node_idx, Vector3s::Zero());
}

Vector3s TwoDScene::getNodePosX(int bucket_idx, int node_idx) const {
  if (m_bucket_activated[bucket_idx]) {
    const scalar dx = getCellSize();
    return m_node_pos[bucket_idx].segment<3>(node_idx * 3) +
           Vector3s(0.0, 0.5, 0.5) * dx;
  } else {
    return nodePosFromBucket(bucket_idx, node_idx, Vector3s(0.0, 0.5, 0.5));
  }
}

Vector3s TwoDScene::getNodePosY(int bucket_idx, int node_idx) const {
  if (m_bucket_activated[bucket_idx]) {
    const scalar dx = getCellSize();
    return m_node_pos[bucket_idx].segment<3>(node_idx * 3) +
           Vector3s(0.5, 0.0, 0.5) * dx;
  } else {
    return nodePosFromBucket(bucket_idx, node_idx, Vector3s(0.5, 0.0, 0.5));
  }
}

Vector3s TwoDScene::getNodePosZ(int bucket_idx, int node_idx) const {
  if (m_bucket_activated[bucket_idx]) {
    const scalar dx = getCellSize();
    return m_node_pos[bucket_idx].segment<3>(node_idx * 3) +
           Vector3s(0.5, 0.5, 0.0) * dx;
  } else {
    return nodePosFromBucket(bucket_idx, node_idx, Vector3s(0.5, 0.5, 0.0));
  }
}

Vector3s TwoDScene::getNodePosP(int bucket_idx, int node_idx) const {
  if (m_bucket_activated[bucket_idx]) {
    const scalar dx = getCellSize();
    return m_node_pos[bucket_idx].segment<3>(node_idx * 3) +
           Vector3s(0.5, 0.5, 0.5) * dx;
  } else {
    return nodePosFromBucket(bucket_idx, node_idx, Vector3s(0.5, 0.5, 0.5));
  }
}

Vector3s TwoDScene::getNodePosEX(int bucket_idx, int node_idx) const {
  if (m_bucket_activated[bucket_idx]) {
    const scalar dx = getCellSize();
    return m_node_pos[bucket_idx].segment<3>(node_idx * 3) +
           Vector3s(0.5, 0.0, 0.0) * dx;
  } else {
    return nodePosFromBucket(bucket_idx, node_idx, Vector3s(0.5, 0.0, 0.0));
  }
}

Vector3s TwoDScene::getNodePosEY(int bucket_idx, int node_idx) const {
  if (m_bucket_activated[bucket_idx]) {
    const scalar dx = getCellSize();
    return m_node_pos[bucket_idx].segment<3>(node_idx * 3) +
           Vector3s(0.0, 0.5, 0.0) * dx;
  } else {
    return nodePosFromBucket(bucket_idx, node_idx, Vector3s(0.0, 0.5, 0.0));
  }
}

Vector3s TwoDScene::getNodePosEZ(int bucket_idx, int node_idx) const {
  if (m_bucket_activated[bucket_idx]) {
    const scalar dx = getCellSize();
    return m_node_pos[bucket_idx].segment<3>(node_idx * 3) +
           Vector3s(0.0, 0.0, 0.5) * dx;
  } else {
    return nodePosFromBucket(bucket_idx, node_idx, Vector3s(0.0, 0.0, 0.5));
  }
}

/*!
 * allocate memory for nodes
 */
void TwoDScene::generateNodes() {
  const scalar dx = getCellSize();

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    // ignore inactivated buckets
    if (!m_bucket_activated[bucket_idx]) return;

    Vector3i bucket_handle = m_particle_buckets.bucket_handle(bucket_idx);

    Vector3s cell_local_corner =
        Vector3s((scalar)bucket_handle(0) * m_bucket_size,
                 (scalar)bucket_handle(1) * m_bucket_size,
                 (scalar)bucket_handle(2) * m_bucket_size) +
        m_grid_mincorner;

    // otherwise generate nodes to fill the bucket
    const int count = m_num_nodes * m_num_nodes * m_num_nodes;
    if ((int)m_node_pos[bucket_idx].size() != count * 3)
      m_node_pos[bucket_idx].resize(count * 3);

    VectorXs& bucket_node_pos = m_node_pos[bucket_idx];

    for (int k = 0; k < m_num_nodes; ++k)
      for (int j = 0; j < m_num_nodes; ++j)
        for (int i = 0; i < m_num_nodes; ++i) {
          int node_idx = k * m_num_nodes * m_num_nodes + j * m_num_nodes + i;
          bucket_node_pos.segment<3>(node_idx * 3) =
              cell_local_corner + Vector3s(i, j, k) * dx;
        }

    VectorXi& bucket_node_idxp_x = m_node_index_pressure_x[bucket_idx];
    VectorXi& bucket_node_idx_solid_phi_x =
        m_node_index_solid_phi_x[bucket_idx];
    bucket_node_idxp_x.resize(count * 4);
    bucket_node_idxp_x.setConstant(-1);
    bucket_node_idx_solid_phi_x.resize(count * 8);
    bucket_node_idx_solid_phi_x.setConstant(-1);

    VectorXi& bucket_node_idxp_y = m_node_index_pressure_y[bucket_idx];
    VectorXi& bucket_node_idx_solid_phi_y =
        m_node_index_solid_phi_y[bucket_idx];
    bucket_node_idxp_y.resize(count * 4);
    bucket_node_idxp_y.setConstant(-1);
    bucket_node_idx_solid_phi_y.resize(count * 8);
    bucket_node_idx_solid_phi_y.setConstant(-1);

    VectorXi& bucket_node_idxp_z = m_node_index_pressure_z[bucket_idx];
    VectorXi& bucket_node_idx_solid_phi_z =
        m_node_index_solid_phi_z[bucket_idx];
    bucket_node_idxp_z.resize(count * 4);
    bucket_node_idxp_z.setConstant(-1);
    bucket_node_idx_solid_phi_z.resize(count * 8);
    bucket_node_idx_solid_phi_z.setConstant(-1);
  });
}

const std::vector<VectorXi>& TwoDScene::getPressureNeighbors() const {
  return m_node_pressure_neighbors;
}

const std::vector<VectorXs>& TwoDScene::getNodeLiquidVolFracCentre() const {
  return m_node_liquid_c_vf;
}

const std::vector<VectorXs>& TwoDScene::getNodeLiquidVolFracU() const {
  return m_node_liquid_u_vf;
}

const std::vector<VectorXs>& TwoDScene::getNodeLiquidVolFracV() const {
  return m_node_liquid_v_vf;
}

const std::vector<VectorXs>& TwoDScene::getNodeLiquidVolFracW() const {
  return m_node_liquid_w_vf;
}

const std::vector<VectorXs>& TwoDScene::getNodeLiquidVolFracEX() const {
  return m_node_liquid_ex_vf;
}

const std::vector<VectorXs>& TwoDScene::getNodeLiquidVolFracEY() const {
  return m_node_liquid_ey_vf;
}

const std::vector<VectorXs>& TwoDScene::getNodeLiquidVolFracEZ() const {
  return m_node_liquid_ez_vf;
}

const std::vector<VectorXs>& TwoDScene::getNodeLiquidPhi() const {
  return m_node_liquid_phi;
}

const std::vector<VectorXi>& TwoDScene::getNodePressureIndexX() const {
  return m_node_index_pressure_x;
}

const std::vector<VectorXi>& TwoDScene::getNodePressureIndexY() const {
  return m_node_index_pressure_y;
}

const std::vector<VectorXi>& TwoDScene::getNodePressureIndexZ() const {
  return m_node_index_pressure_z;
}

const std::vector<VectorXs>& TwoDScene::getNodeSolidWeightX() const {
  return m_node_solid_weight_x;
}

const std::vector<VectorXs>& TwoDScene::getNodeSolidWeightY() const {
  return m_node_solid_weight_y;
}

const std::vector<VectorXs>& TwoDScene::getNodeSolidWeightZ() const {
  return m_node_solid_weight_z;
}

const std::vector<VectorXs>& TwoDScene::getNodeSolidVelX() const {
  return m_node_solid_vel_x;
}

const std::vector<VectorXs>& TwoDScene::getNodeSolidVelY() const {
  return m_node_solid_vel_y;
}

const std::vector<VectorXs>& TwoDScene::getNodeSolidVelZ() const {
  return m_node_solid_vel_z;
}

void TwoDScene::setSimInfo(const SimInfo& info) { m_sim_info = info; }

const std::vector<VectorXs>& TwoDScene::getNodePorePressureP() const {
  return m_node_pore_pressure_p;
}

std::vector<VectorXs>& TwoDScene::getNodePorePressureP() {
  return m_node_pore_pressure_p;
}

const std::vector<VectorXs>& TwoDScene::getNodeSaturationP() const {
  return m_node_sat_p;
}

std::vector<VectorXs>& TwoDScene::getNodeSaturationP() { return m_node_sat_p; }

const std::vector<VectorXs>& TwoDScene::getNodeSaturationX() const {
  return m_node_sat_x;
}

std::vector<VectorXs>& TwoDScene::getNodeSaturationX() { return m_node_sat_x; }

const std::vector<VectorXs>& TwoDScene::getNodeSaturationY() const {
  return m_node_sat_y;
}

std::vector<VectorXs>& TwoDScene::getNodeSaturationY() { return m_node_sat_y; }

const std::vector<VectorXs>& TwoDScene::getNodeSaturationZ() const {
  return m_node_sat_z;
}

std::vector<VectorXs>& TwoDScene::getNodeSaturationZ() { return m_node_sat_z; }

const std::vector<VectorXs>& TwoDScene::getNodePsiX() const {
  return m_node_psi_x;
}

std::vector<VectorXs>& TwoDScene::getNodePsiX() { return m_node_psi_x; }

const std::vector<VectorXs>& TwoDScene::getNodePsiY() const {
  return m_node_psi_y;
}

std::vector<VectorXs>& TwoDScene::getNodePsiY() { return m_node_psi_y; }

const std::vector<VectorXs>& TwoDScene::getNodePsiZ() const {
  return m_node_psi_z;
}

std::vector<VectorXs>& TwoDScene::getNodePsiZ() { return m_node_psi_z; }

scalar TwoDScene::getMaxVelocity() const {
  scalar max_vel = 0.0;
  const int num_elasto = getNumSoftElastoParticles();
  for (int i = 0; i < num_elasto; ++i) {
    max_vel = std::max(max_vel, m_v.segment<3>(i * 4).squaredNorm());
  }
  return sqrt(max_vel);
}

scalar TwoDScene::getMaxFluidVelocity() const {
  scalar max_vel = 0.0;
  const int num_fluid = getNumFluidParticles();
  for (int i = 0; i < num_fluid; ++i) {
    max_vel =
        std::max(max_vel, m_fluid_v.segment<3>(m_fluids[i] * 4).squaredNorm());
  }
  return sqrt(max_vel);
}

const std::vector<VectorXi>& TwoDScene::getNodeIndexEdgeX() const {
  return m_node_index_edge_x;
}

const std::vector<VectorXi>& TwoDScene::getNodeIndexEdgeY() const {
  return m_node_index_edge_y;
}

const std::vector<VectorXi>& TwoDScene::getNodeIndexEdgeZ() const {
  return m_node_index_edge_z;
}

/*!
 * connect edges to neighbor nodes
 */
void TwoDScene::connectEdgeNodes() {
  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    Vector3i bucket_handle = m_particle_buckets.bucket_handle(bucket_idx);

    VectorXi& bucket_node_idx_ex = m_node_index_edge_x[bucket_idx];
    VectorXi& bucket_node_idx_ey = m_node_index_edge_y[bucket_idx];
    VectorXi& bucket_node_idx_ez = m_node_index_edge_z[bucket_idx];

    if (!m_bucket_activated[bucket_idx]) return;

    for (int k = 0; k < m_num_nodes; ++k)
      for (int j = 0; j < m_num_nodes; ++j)
        for (int i = 0; i < m_num_nodes; ++i) {
          int node_idx = k * m_num_nodes * m_num_nodes + j * m_num_nodes + i;

          // back, front (edge-y)
          for (int r = 0; r < 2; ++r) {
            Vector3i node_bucket_handle = bucket_handle;
            Vector3i ey_local_idx = Vector3i(i, j, k + r);

            if (ey_local_idx(2) >= m_num_nodes) {
              node_bucket_handle(2)++;
              ey_local_idx(2) -= m_num_nodes;
            }

            if (node_bucket_handle(0) < 0 ||
                node_bucket_handle(0) >= m_particle_buckets.dim_size(0) ||
                node_bucket_handle(1) < 0 ||
                node_bucket_handle(1) >= m_particle_buckets.dim_size(1) ||
                node_bucket_handle(2) < 0 ||
                node_bucket_handle(2) >= m_particle_buckets.dim_size(2)) {
              bucket_node_idx_ex(node_idx * 8 + r * 2 + 0) = -1;
              bucket_node_idx_ex(node_idx * 8 + r * 2 + 1) = -1;
            } else {
              const int ey_idx = ey_local_idx(2) * m_num_nodes * m_num_nodes +
                                 ey_local_idx(1) * m_num_nodes +
                                 ey_local_idx(0);
              int nb_bucket_idx =
                  m_particle_buckets.bucket_index(node_bucket_handle);
              if (!m_bucket_activated[nb_bucket_idx]) {
                bucket_node_idx_ex(node_idx * 8 + r * 2 + 0) = -1;
                bucket_node_idx_ex(node_idx * 8 + r * 2 + 1) = -1;
              } else {
                bucket_node_idx_ex(node_idx * 8 + r * 2 + 0) = nb_bucket_idx;
                bucket_node_idx_ex(node_idx * 8 + r * 2 + 1) = ey_idx;
              }
            }
          }

          // bottom, top (edge-z)
          for (int r = 0; r < 2; ++r) {
            Vector3i node_bucket_handle = bucket_handle;
            Vector3i ez_local_idx = Vector3i(i, j + r, k);

            if (ez_local_idx(1) >= m_num_nodes) {
              node_bucket_handle(1)++;
              ez_local_idx(1) -= m_num_nodes;
            }

            if (node_bucket_handle(0) < 0 ||
                node_bucket_handle(0) >= m_particle_buckets.dim_size(0) ||
                node_bucket_handle(1) < 0 ||
                node_bucket_handle(1) >= m_particle_buckets.dim_size(1) ||
                node_bucket_handle(2) < 0 ||
                node_bucket_handle(2) >= m_particle_buckets.dim_size(2)) {
              bucket_node_idx_ex(node_idx * 8 + 4 + r * 2 + 0) = -1;
              bucket_node_idx_ex(node_idx * 8 + 4 + r * 2 + 1) = -1;
            } else {
              const int ez_idx = ez_local_idx(2) * m_num_nodes * m_num_nodes +
                                 ez_local_idx(1) * m_num_nodes +
                                 ez_local_idx(0);
              int nb_bucket_idx =
                  m_particle_buckets.bucket_index(node_bucket_handle);
              if (!m_bucket_activated[nb_bucket_idx]) {
                bucket_node_idx_ex(node_idx * 8 + 4 + r * 2 + 0) = -1;
                bucket_node_idx_ex(node_idx * 8 + 4 + r * 2 + 1) = -1;
              } else {
                bucket_node_idx_ex(node_idx * 8 + 4 + r * 2 + 0) =
                    nb_bucket_idx;
                bucket_node_idx_ex(node_idx * 8 + 4 + r * 2 + 1) = ez_idx;
              }
            }
          }

          // back, front (edge-x)
          for (int r = 0; r < 2; ++r) {
            Vector3i node_bucket_handle = bucket_handle;
            Vector3i ex_local_idx = Vector3i(i, j, k + r);

            if (ex_local_idx(2) >= m_num_nodes) {
              node_bucket_handle(2)++;
              ex_local_idx(2) -= m_num_nodes;
            }

            if (node_bucket_handle(0) < 0 ||
                node_bucket_handle(0) >= m_particle_buckets.dim_size(0) ||
                node_bucket_handle(1) < 0 ||
                node_bucket_handle(1) >= m_particle_buckets.dim_size(1) ||
                node_bucket_handle(2) < 0 ||
                node_bucket_handle(2) >= m_particle_buckets.dim_size(2)) {
              bucket_node_idx_ey(node_idx * 8 + r * 2 + 0) = -1;
              bucket_node_idx_ey(node_idx * 8 + r * 2 + 1) = -1;
            } else {
              const int ex_idx = ex_local_idx(2) * m_num_nodes * m_num_nodes +
                                 ex_local_idx(1) * m_num_nodes +
                                 ex_local_idx(0);
              int nb_bucket_idx =
                  m_particle_buckets.bucket_index(node_bucket_handle);
              if (!m_bucket_activated[nb_bucket_idx]) {
                bucket_node_idx_ey(node_idx * 8 + r * 2 + 0) = -1;
                bucket_node_idx_ey(node_idx * 8 + r * 2 + 1) = -1;
              } else {
                bucket_node_idx_ey(node_idx * 8 + r * 2 + 0) = nb_bucket_idx;
                bucket_node_idx_ey(node_idx * 8 + r * 2 + 1) = ex_idx;
              }
            }
          }

          // left, right (edge-z)
          for (int r = 0; r < 2; ++r) {
            Vector3i node_bucket_handle = bucket_handle;
            Vector3i ez_local_idx = Vector3i(i + r, j, k);

            if (ez_local_idx(0) >= m_num_nodes) {
              node_bucket_handle(0)++;
              ez_local_idx(0) -= m_num_nodes;
            }

            if (node_bucket_handle(0) < 0 ||
                node_bucket_handle(0) >= m_particle_buckets.dim_size(0) ||
                node_bucket_handle(1) < 0 ||
                node_bucket_handle(1) >= m_particle_buckets.dim_size(1) ||
                node_bucket_handle(2) < 0 ||
                node_bucket_handle(2) >= m_particle_buckets.dim_size(2)) {
              bucket_node_idx_ey(node_idx * 8 + 4 + r * 2 + 0) = -1;
              bucket_node_idx_ey(node_idx * 8 + 4 + r * 2 + 1) = -1;
            } else {
              const int ez_idx = ez_local_idx(2) * m_num_nodes * m_num_nodes +
                                 ez_local_idx(1) * m_num_nodes +
                                 ez_local_idx(0);
              int nb_bucket_idx =
                  m_particle_buckets.bucket_index(node_bucket_handle);

              if (!m_bucket_activated[nb_bucket_idx]) {
                bucket_node_idx_ey(node_idx * 8 + 4 + r * 2 + 0) = -1;
                bucket_node_idx_ey(node_idx * 8 + 4 + r * 2 + 1) = -1;
              } else {
                bucket_node_idx_ey(node_idx * 8 + 4 + r * 2 + 0) =
                    nb_bucket_idx;
                bucket_node_idx_ey(node_idx * 8 + 4 + r * 2 + 1) = ez_idx;
              }
            }
          }

          // bottom, top (edge-x)
          for (int r = 0; r < 2; ++r) {
            Vector3i node_bucket_handle = bucket_handle;
            Vector3i ex_local_idx = Vector3i(i, j + r, k);

            if (ex_local_idx(1) >= m_num_nodes) {
              node_bucket_handle(1)++;
              ex_local_idx(1) -= m_num_nodes;
            }

            if (node_bucket_handle(0) < 0 ||
                node_bucket_handle(0) >= m_particle_buckets.dim_size(0) ||
                node_bucket_handle(1) < 0 ||
                node_bucket_handle(1) >= m_particle_buckets.dim_size(1) ||
                node_bucket_handle(2) < 0 ||
                node_bucket_handle(2) >= m_particle_buckets.dim_size(2)) {
              bucket_node_idx_ez(node_idx * 8 + r * 2 + 0) = -1;
              bucket_node_idx_ez(node_idx * 8 + r * 2 + 1) = -1;
            } else {
              const int ex_idx = ex_local_idx(2) * m_num_nodes * m_num_nodes +
                                 ex_local_idx(1) * m_num_nodes +
                                 ex_local_idx(0);
              int nb_bucket_idx =
                  m_particle_buckets.bucket_index(node_bucket_handle);

              if (!m_bucket_activated[nb_bucket_idx]) {
                bucket_node_idx_ez(node_idx * 8 + r * 2 + 0) = -1;
                bucket_node_idx_ez(node_idx * 8 + r * 2 + 1) = -1;
              } else {
                bucket_node_idx_ez(node_idx * 8 + r * 2 + 0) = nb_bucket_idx;
                bucket_node_idx_ez(node_idx * 8 + r * 2 + 1) = ex_idx;
              }
            }
          }

          // left, right (edge-y)
          for (int r = 0; r < 2; ++r) {
            Vector3i node_bucket_handle = bucket_handle;
            Vector3i ey_local_idx = Vector3i(i + r, j, k);

            if (ey_local_idx(0) >= m_num_nodes) {
              node_bucket_handle(0)++;
              ey_local_idx(0) -= m_num_nodes;
            }

            if (node_bucket_handle(0) < 0 ||
                node_bucket_handle(0) >= m_particle_buckets.dim_size(0) ||
                node_bucket_handle(1) < 0 ||
                node_bucket_handle(1) >= m_particle_buckets.dim_size(1) ||
                node_bucket_handle(2) < 0 ||
                node_bucket_handle(2) >= m_particle_buckets.dim_size(2)) {
              bucket_node_idx_ez(node_idx * 8 + 4 + r * 2 + 0) = -1;
              bucket_node_idx_ez(node_idx * 8 + 4 + r * 2 + 1) = -1;
            } else {
              const int ey_idx = ey_local_idx(2) * m_num_nodes * m_num_nodes +
                                 ey_local_idx(1) * m_num_nodes +
                                 ey_local_idx(0);
              int nb_bucket_idx =
                  m_particle_buckets.bucket_index(node_bucket_handle);

              if (!m_bucket_activated[nb_bucket_idx]) {
                bucket_node_idx_ez(node_idx * 8 + 4 + r * 2 + 0) = -1;
                bucket_node_idx_ez(node_idx * 8 + 4 + r * 2 + 1) = -1;
              } else {
                bucket_node_idx_ez(node_idx * 8 + 4 + r * 2 + 0) =
                    nb_bucket_idx;
                bucket_node_idx_ez(node_idx * 8 + 4 + r * 2 + 1) = ey_idx;
              }
            }
          }
        }
  });
}

void TwoDScene::connectSolidPhiNodes() {
  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    Vector3i bucket_handle = m_particle_buckets.bucket_handle(bucket_idx);

    VectorXi& bucket_node_idx_solid_phi_x =
        m_node_index_solid_phi_x[bucket_idx];
    VectorXi& bucket_node_idx_solid_phi_y =
        m_node_index_solid_phi_y[bucket_idx];
    VectorXi& bucket_node_idx_solid_phi_z =
        m_node_index_solid_phi_z[bucket_idx];

    if (!m_bucket_activated[bucket_idx]) return;

    for (int k = 0; k < m_num_nodes; ++k)
      for (int j = 0; j < m_num_nodes; ++j)
        for (int i = 0; i < m_num_nodes; ++i) {
          int node_idx = k * m_num_nodes * m_num_nodes + j * m_num_nodes + i;

          for (int r = 0; r < 2; ++r)
            for (int s = 0; s < 2; ++s) {
              Vector3i node_bucket_handle_x = bucket_handle;

              Vector3i sphi_local_idx = Vector3i(i, j + s, k + r);

              if (sphi_local_idx(1) >= m_num_nodes) {
                node_bucket_handle_x(1)++;
                sphi_local_idx(1) -= m_num_nodes;
              }

              if (sphi_local_idx(2) >= m_num_nodes) {
                node_bucket_handle_x(2)++;
                sphi_local_idx(2) -= m_num_nodes;
              }

              if (node_bucket_handle_x(0) < 0 ||
                  node_bucket_handle_x(0) >= m_particle_buckets.dim_size(0) ||
                  node_bucket_handle_x(1) < 0 ||
                  node_bucket_handle_x(1) >= m_particle_buckets.dim_size(1) ||
                  node_bucket_handle_x(2) < 0 ||
                  node_bucket_handle_x(2) >= m_particle_buckets.dim_size(2)) {
                bucket_node_idx_solid_phi_x(node_idx * 8 + (r * 2 + s) * 2 +
                                            0) = -1;
                bucket_node_idx_solid_phi_x(node_idx * 8 + (r * 2 + s) * 2 +
                                            1) = -1;
              } else {
                const int sphi_idx =
                    sphi_local_idx(2) * m_num_nodes * m_num_nodes +
                    sphi_local_idx(1) * m_num_nodes + sphi_local_idx(0);

                int nb_bucket_idx =
                    m_particle_buckets.bucket_index(node_bucket_handle_x);
                if (!m_bucket_activated[nb_bucket_idx]) {
                  bucket_node_idx_solid_phi_x(node_idx * 8 + (r * 2 + s) * 2 +
                                              0) = -1;
                  bucket_node_idx_solid_phi_x(node_idx * 8 + (r * 2 + s) * 2 +
                                              1) = -1;
                } else {
                  bucket_node_idx_solid_phi_x(node_idx * 8 + (r * 2 + s) * 2 +
                                              0) = nb_bucket_idx;
                  bucket_node_idx_solid_phi_x(node_idx * 8 + (r * 2 + s) * 2 +
                                              1) = sphi_idx;
                }
              }
            }

          for (int r = 0; r < 2; ++r)
            for (int s = 0; s < 2; ++s) {
              Vector3i node_bucket_handle_y = bucket_handle;

              Vector3i sphi_local_idx = Vector3i(i + r, j, k + s);

              if (sphi_local_idx(0) >= m_num_nodes) {
                node_bucket_handle_y(0)++;
                sphi_local_idx(0) -= m_num_nodes;
              }

              if (sphi_local_idx(2) >= m_num_nodes) {
                node_bucket_handle_y(2)++;
                sphi_local_idx(2) -= m_num_nodes;
              }

              if (node_bucket_handle_y(0) < 0 ||
                  node_bucket_handle_y(0) >= m_particle_buckets.dim_size(0) ||
                  node_bucket_handle_y(1) < 0 ||
                  node_bucket_handle_y(1) >= m_particle_buckets.dim_size(1) ||
                  node_bucket_handle_y(2) < 0 ||
                  node_bucket_handle_y(2) >= m_particle_buckets.dim_size(2)) {
                bucket_node_idx_solid_phi_y(node_idx * 8 + (r * 2 + s) * 2 +
                                            0) = -1;
                bucket_node_idx_solid_phi_y(node_idx * 8 + (r * 2 + s) * 2 +
                                            1) = -1;
              } else {
                const int sphi_idx =
                    sphi_local_idx(2) * m_num_nodes * m_num_nodes +
                    sphi_local_idx(1) * m_num_nodes + sphi_local_idx(0);

                int nb_bucket_idx =
                    m_particle_buckets.bucket_index(node_bucket_handle_y);
                if (!m_bucket_activated[nb_bucket_idx]) {
                  bucket_node_idx_solid_phi_y(node_idx * 8 + (r * 2 + s) * 2 +
                                              0) = -1;
                  bucket_node_idx_solid_phi_y(node_idx * 8 + (r * 2 + s) * 2 +
                                              1) = -1;
                } else {
                  bucket_node_idx_solid_phi_y(node_idx * 8 + (r * 2 + s) * 2 +
                                              0) = nb_bucket_idx;
                  bucket_node_idx_solid_phi_y(node_idx * 8 + (r * 2 + s) * 2 +
                                              1) = sphi_idx;
                }
              }
            }

          for (int r = 0; r < 2; ++r)
            for (int s = 0; s < 2; ++s) {
              Vector3i node_bucket_handle_z = bucket_handle;

              Vector3i sphi_local_idx = Vector3i(i + r, j + s, k);

              if (sphi_local_idx(0) >= m_num_nodes) {
                node_bucket_handle_z(0)++;
                sphi_local_idx(0) -= m_num_nodes;
              }

              if (sphi_local_idx(1) >= m_num_nodes) {
                node_bucket_handle_z(1)++;
                sphi_local_idx(1) -= m_num_nodes;
              }

              if (node_bucket_handle_z(0) < 0 ||
                  node_bucket_handle_z(0) >= m_particle_buckets.dim_size(0) ||
                  node_bucket_handle_z(1) < 0 ||
                  node_bucket_handle_z(1) >= m_particle_buckets.dim_size(1) ||
                  node_bucket_handle_z(2) < 0 ||
                  node_bucket_handle_z(2) >= m_particle_buckets.dim_size(2)) {
                bucket_node_idx_solid_phi_z(node_idx * 8 + (r * 2 + s) * 2 +
                                            0) = -1;
                bucket_node_idx_solid_phi_z(node_idx * 8 + (r * 2 + s) * 2 +
                                            1) = -1;
              } else {
                const int sphi_idx =
                    sphi_local_idx(2) * m_num_nodes * m_num_nodes +
                    sphi_local_idx(1) * m_num_nodes + sphi_local_idx(0);

                int nb_bucket_idx =
                    m_particle_buckets.bucket_index(node_bucket_handle_z);
                if (!m_bucket_activated[nb_bucket_idx]) {
                  bucket_node_idx_solid_phi_z(node_idx * 8 + (r * 2 + s) * 2 +
                                              0) = -1;
                  bucket_node_idx_solid_phi_z(node_idx * 8 + (r * 2 + s) * 2 +
                                              1) = -1;
                } else {
                  bucket_node_idx_solid_phi_z(node_idx * 8 + (r * 2 + s) * 2 +
                                              0) = nb_bucket_idx;
                  bucket_node_idx_solid_phi_z(node_idx * 8 + (r * 2 + s) * 2 +
                                              1) = sphi_idx;
                }
              }
            }
        }
  });
}

/*!
 * connect pressure nodes with X- Y- and Z- nodes
 */
void TwoDScene::connectPressureNodes() {
  const Vector3i ppdir[] = {
      Vector3i(-1, 0, 0),  Vector3i(1, 0, 0),   Vector3i(0, -1, 0),
      Vector3i(0, 1, 0),   Vector3i(0, 0, -1),  Vector3i(0, 0, 1),
      Vector3i(-1, -1, 0), Vector3i(1, -1, 0),  Vector3i(-1, 1, 0),
      Vector3i(1, 1, 0),   Vector3i(-1, 0, -1), Vector3i(1, 0, -1),
      Vector3i(0, -1, -1), Vector3i(0, 1, -1),  Vector3i(-1, 0, 1),
      Vector3i(1, 0, 1),   Vector3i(0, -1, 1),  Vector3i(0, 1, 1)};

  m_particle_buckets.for_each_bucket_colored([&](int bucket_idx) {
    Vector3i bucket_handle = m_particle_buckets.bucket_handle(bucket_idx);
    if (!m_bucket_activated[bucket_idx]) return;

    VectorXi& bucket_node_pressure_neighbors =
        m_node_pressure_neighbors[bucket_idx];
    VectorXi& bucket_node_pp_neighbors = m_node_pp_neighbors[bucket_idx];

    const int count_p = getNumNodes(bucket_idx);
    bucket_node_pressure_neighbors.resize(count_p * 6 * 2);
    bucket_node_pp_neighbors.resize(count_p * 18 * 2);

    bucket_node_pressure_neighbors.setConstant(-1);
    bucket_node_pp_neighbors.setConstant(-1);

    for (int k = 0; k < m_num_nodes; ++k)
      for (int j = 0; j < m_num_nodes; ++j)
        for (int i = 0; i < m_num_nodes; ++i) {
          int node_idx = k * m_num_nodes * m_num_nodes + j * m_num_nodes + i;

          for (int r = 0; r < 18; ++r) {
            Vector3i node_bucket_handle_p = bucket_handle;
            Vector3i mac_local_idx = Vector3i(i, j, k) + ppdir[r];

            for (int s = 0; s < 3; ++s) {
              if (mac_local_idx(s) >= m_num_nodes) {
                node_bucket_handle_p(s)++;
                mac_local_idx(s) -= m_num_nodes;
              } else if (mac_local_idx(s) < 0) {
                node_bucket_handle_p(s)--;
                mac_local_idx(s) += m_num_nodes;
              }
            }

            if (!m_particle_buckets.has_bucket(node_bucket_handle_p)) continue;

            const int nb_bucket_idx =
                m_particle_buckets.bucket_index(node_bucket_handle_p);

            if (!m_bucket_activated[nb_bucket_idx]) continue;

            const int mac_idx = mac_local_idx(2) * m_num_nodes * m_num_nodes +
                                mac_local_idx(1) * m_num_nodes +
                                mac_local_idx(0);

            bucket_node_pp_neighbors[node_idx * 36 + r * 2 + 0] = nb_bucket_idx;
            bucket_node_pp_neighbors[node_idx * 36 + r * 2 + 1] = mac_idx;
          }

          for (int r = 0; r < 2; ++r) {
            Vector3i node_bucket_handle_x = bucket_handle;

            Vector3i mac_local_idx = Vector3i(i + r, j, k);

            if (mac_local_idx(0) >= m_num_nodes) {
              node_bucket_handle_x(0)++;
              mac_local_idx(0) -= m_num_nodes;
            }

            if (!m_particle_buckets.has_bucket(node_bucket_handle_x)) continue;

            int nb_bucket_idx =
                m_particle_buckets.bucket_index(node_bucket_handle_x);

            if (!m_bucket_activated[nb_bucket_idx]) continue;

            const int mac_idx = mac_local_idx(2) * m_num_nodes * m_num_nodes +
                                mac_local_idx(1) * m_num_nodes +
                                mac_local_idx(0);

            bucket_node_pressure_neighbors[node_idx * 12 + r * 2 + 0] =
                nb_bucket_idx;
            bucket_node_pressure_neighbors[node_idx * 12 + r * 2 + 1] = mac_idx;

            VectorXi& nb_node_idxp = m_node_index_pressure_x[nb_bucket_idx];
            nb_node_idxp[mac_idx * 4 + (1 - r) * 2 + 0] = bucket_idx;
            nb_node_idxp[mac_idx * 4 + (1 - r) * 2 + 1] = node_idx;
          }

          for (int r = 0; r < 2; ++r) {
            Vector3i node_bucket_handle_y = bucket_handle;

            Vector3i mac_local_idx = Vector3i(i, j + r, k);

            if (mac_local_idx(1) >= m_num_nodes) {
              node_bucket_handle_y(1)++;
              mac_local_idx(1) -= m_num_nodes;
            }

            if (!m_particle_buckets.has_bucket(node_bucket_handle_y)) continue;

            int nb_bucket_idx =
                m_particle_buckets.bucket_index(node_bucket_handle_y);

            if (!m_bucket_activated[nb_bucket_idx]) continue;

            int mac_idx = mac_local_idx(2) * m_num_nodes * m_num_nodes +
                          mac_local_idx(1) * m_num_nodes + mac_local_idx(0);

            bucket_node_pressure_neighbors[node_idx * 12 + 4 + r * 2 + 0] =
                nb_bucket_idx;
            bucket_node_pressure_neighbors[node_idx * 12 + 4 + r * 2 + 1] =
                mac_idx;

            VectorXi& nb_node_idxp = m_node_index_pressure_y[nb_bucket_idx];
            nb_node_idxp[mac_idx * 4 + (1 - r) * 2 + 0] = bucket_idx;
            nb_node_idxp[mac_idx * 4 + (1 - r) * 2 + 1] = node_idx;
          }

          for (int r = 0; r < 2; ++r) {
            Vector3i node_bucket_handle_z = bucket_handle;

            Vector3i mac_local_idx = Vector3i(i, j, k + r);

            if (mac_local_idx(2) >= m_num_nodes) {
              node_bucket_handle_z(2)++;
              mac_local_idx(2) -= m_num_nodes;
            }

            if (!m_particle_buckets.has_bucket(node_bucket_handle_z)) continue;

            int nb_bucket_idx =
                m_particle_buckets.bucket_index(node_bucket_handle_z);

            if (!m_bucket_activated[nb_bucket_idx]) continue;

            int mac_idx = mac_local_idx(2) * m_num_nodes * m_num_nodes +
                          mac_local_idx(1) * m_num_nodes + mac_local_idx(0);

            bucket_node_pressure_neighbors[node_idx * 12 + 8 + r * 2 + 0] =
                nb_bucket_idx;
            bucket_node_pressure_neighbors[node_idx * 12 + 8 + r * 2 + 1] =
                mac_idx;

            VectorXi& nb_node_idxp = m_node_index_pressure_z[nb_bucket_idx];
            nb_node_idxp[mac_idx * 4 + (1 - r) * 2 + 0] = bucket_idx;
            nb_node_idxp[mac_idx * 4 + (1 - r) * 2 + 1] = node_idx;
          }
        }
  });
}

/*!
 * allocate attributes to be stored on nodes
 */
void TwoDScene::postAllocateNodes() {
  const int num_buckets = m_particle_buckets.size();
  // check allocation
  if ((int)m_node_mass_x.size() != num_buckets)
    m_node_mass_x.resize(num_buckets);
  if ((int)m_node_sat_x.size() != num_buckets) m_node_sat_x.resize(num_buckets);
  if ((int)m_node_psi_x.size() != num_buckets) m_node_psi_x.resize(num_buckets);
  if ((int)m_node_vel_x.size() != num_buckets) m_node_vel_x.resize(num_buckets);
  if ((int)m_node_vol_x.size() != num_buckets) m_node_vol_x.resize(num_buckets);
  if ((int)m_node_shape_factor_x.size() != num_buckets)
    m_node_shape_factor_x.resize(num_buckets);
  if ((int)m_node_raw_weight_x.size() != num_buckets)
    m_node_raw_weight_x.resize(num_buckets);
  if ((int)m_node_orientation_x.size() != num_buckets)
    m_node_orientation_x.resize(num_buckets);

  if ((int)m_node_mass_y.size() != num_buckets)
    m_node_mass_y.resize(num_buckets);
  if ((int)m_node_sat_y.size() != num_buckets) m_node_sat_y.resize(num_buckets);
  if ((int)m_node_psi_y.size() != num_buckets) m_node_psi_y.resize(num_buckets);
  if ((int)m_node_vel_y.size() != num_buckets) m_node_vel_y.resize(num_buckets);
  if ((int)m_node_vol_y.size() != num_buckets) m_node_vol_y.resize(num_buckets);
  if ((int)m_node_shape_factor_y.size() != num_buckets)
    m_node_shape_factor_y.resize(num_buckets);
  if ((int)m_node_raw_weight_y.size() != num_buckets)
    m_node_raw_weight_y.resize(num_buckets);
  if ((int)m_node_orientation_y.size() != num_buckets)
    m_node_orientation_y.resize(num_buckets);

  if ((int)m_node_mass_z.size() != num_buckets)
    m_node_mass_z.resize(num_buckets);
  if ((int)m_node_sat_z.size() != num_buckets) m_node_sat_z.resize(num_buckets);
  if ((int)m_node_psi_z.size() != num_buckets) m_node_psi_z.resize(num_buckets);
  if ((int)m_node_vel_z.size() != num_buckets) m_node_vel_z.resize(num_buckets);
  if ((int)m_node_vol_z.size() != num_buckets) m_node_vol_z.resize(num_buckets);
  if ((int)m_node_shape_factor_z.size() != num_buckets)
    m_node_shape_factor_z.resize(num_buckets);
  if ((int)m_node_raw_weight_z.size() != num_buckets)
    m_node_raw_weight_z.resize(num_buckets);
  if ((int)m_node_orientation_z.size() != num_buckets)
    m_node_orientation_z.resize(num_buckets);

  if ((int)m_node_mass_fluid_x.size() != num_buckets)
    m_node_mass_fluid_x.resize(num_buckets);
  if ((int)m_node_vel_fluid_x.size() != num_buckets)
    m_node_vel_fluid_x.resize(num_buckets);
  if ((int)m_node_vol_fluid_x.size() != num_buckets)
    m_node_vol_fluid_x.resize(num_buckets);
  if ((int)m_node_vol_pure_fluid_x.size() != num_buckets)
    m_node_vol_pure_fluid_x.resize(num_buckets);

  if ((int)m_node_mass_fluid_y.size() != num_buckets)
    m_node_mass_fluid_y.resize(num_buckets);
  if ((int)m_node_vel_fluid_y.size() != num_buckets)
    m_node_vel_fluid_y.resize(num_buckets);
  if ((int)m_node_vol_fluid_y.size() != num_buckets)
    m_node_vol_fluid_y.resize(num_buckets);
  if ((int)m_node_vol_pure_fluid_y.size() != num_buckets)
    m_node_vol_pure_fluid_y.resize(num_buckets);

  if ((int)m_node_mass_fluid_z.size() != num_buckets)
    m_node_mass_fluid_z.resize(num_buckets);
  if ((int)m_node_vel_fluid_z.size() != num_buckets)
    m_node_vel_fluid_z.resize(num_buckets);
  if ((int)m_node_vol_fluid_z.size() != num_buckets)
    m_node_vol_fluid_z.resize(num_buckets);
  if ((int)m_node_vol_pure_fluid_z.size() != num_buckets)
    m_node_vol_pure_fluid_z.resize(num_buckets);

  if ((int)m_node_solid_phi.size() != num_buckets)
    m_node_solid_phi.resize(num_buckets);

  if ((int)m_node_solid_vel_x.size() != num_buckets)
    m_node_solid_vel_x.resize(num_buckets);
  if ((int)m_node_solid_vel_y.size() != num_buckets)
    m_node_solid_vel_y.resize(num_buckets);
  if ((int)m_node_solid_vel_z.size() != num_buckets)
    m_node_solid_vel_z.resize(num_buckets);

  if ((int)m_node_liquid_valid_x.size() != num_buckets)
    m_node_liquid_valid_x.resize(num_buckets);
  if ((int)m_node_liquid_valid_y.size() != num_buckets)
    m_node_liquid_valid_y.resize(num_buckets);
  if ((int)m_node_liquid_valid_z.size() != num_buckets)
    m_node_liquid_valid_z.resize(num_buckets);

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    const int num_nodes = getNumNodes(bucket_idx);

    if (m_node_mass_x[bucket_idx].size() != num_nodes)
      m_node_mass_x[bucket_idx].resize(num_nodes);
    if (m_node_vel_x[bucket_idx].size() != num_nodes)
      m_node_vel_x[bucket_idx].resize(num_nodes);
    if (m_node_vol_x[bucket_idx].size() != num_nodes)
      m_node_vol_x[bucket_idx].resize(num_nodes);
    if (m_node_sat_x[bucket_idx].size() != num_nodes)
      m_node_sat_x[bucket_idx].resize(num_nodes);
    if (m_node_psi_x[bucket_idx].size() != num_nodes)
      m_node_psi_x[bucket_idx].resize(num_nodes);
    if (m_node_shape_factor_x[bucket_idx].size() != num_nodes)
      m_node_shape_factor_x[bucket_idx].resize(num_nodes);
    if (m_node_raw_weight_x[bucket_idx].size() != num_nodes)
      m_node_raw_weight_x[bucket_idx].resize(num_nodes);
    if (m_node_orientation_x[bucket_idx].size() != num_nodes * 3)
      m_node_orientation_x[bucket_idx].resize(num_nodes * 3);

    if (m_node_mass_y[bucket_idx].size() != num_nodes)
      m_node_mass_y[bucket_idx].resize(num_nodes);
    if (m_node_vel_y[bucket_idx].size() != num_nodes)
      m_node_vel_y[bucket_idx].resize(num_nodes);
    if (m_node_vol_y[bucket_idx].size() != num_nodes)
      m_node_vol_y[bucket_idx].resize(num_nodes);
    if (m_node_sat_y[bucket_idx].size() != num_nodes)
      m_node_sat_y[bucket_idx].resize(num_nodes);
    if (m_node_psi_y[bucket_idx].size() != num_nodes)
      m_node_psi_y[bucket_idx].resize(num_nodes);
    if (m_node_shape_factor_y[bucket_idx].size() != num_nodes)
      m_node_shape_factor_y[bucket_idx].resize(num_nodes);
    if (m_node_raw_weight_y[bucket_idx].size() != num_nodes)
      m_node_raw_weight_y[bucket_idx].resize(num_nodes);
    if (m_node_orientation_y[bucket_idx].size() != num_nodes * 3)
      m_node_orientation_y[bucket_idx].resize(num_nodes * 3);

    if (m_node_mass_z[bucket_idx].size() != num_nodes)
      m_node_mass_z[bucket_idx].resize(num_nodes);
    if (m_node_vel_z[bucket_idx].size() != num_nodes)
      m_node_vel_z[bucket_idx].resize(num_nodes);
    if (m_node_vol_z[bucket_idx].size() != num_nodes)
      m_node_vol_z[bucket_idx].resize(num_nodes);
    if (m_node_sat_z[bucket_idx].size() != num_nodes)
      m_node_sat_z[bucket_idx].resize(num_nodes);
    if (m_node_psi_z[bucket_idx].size() != num_nodes)
      m_node_psi_z[bucket_idx].resize(num_nodes);
    if (m_node_shape_factor_z[bucket_idx].size() != num_nodes)
      m_node_shape_factor_z[bucket_idx].resize(num_nodes);
    if (m_node_raw_weight_z[bucket_idx].size() != num_nodes)
      m_node_raw_weight_z[bucket_idx].resize(num_nodes);
    if (m_node_orientation_z[bucket_idx].size() != num_nodes * 3)
      m_node_orientation_z[bucket_idx].resize(num_nodes * 3);

    if (m_node_mass_fluid_x[bucket_idx].size() != num_nodes)
      m_node_mass_fluid_x[bucket_idx].resize(num_nodes);
    if (m_node_vel_fluid_x[bucket_idx].size() != num_nodes)
      m_node_vel_fluid_x[bucket_idx].resize(num_nodes);
    if (m_node_vol_fluid_x[bucket_idx].size() != num_nodes)
      m_node_vol_fluid_x[bucket_idx].resize(num_nodes);
    if (m_node_vol_pure_fluid_x[bucket_idx].size() != num_nodes)
      m_node_vol_pure_fluid_x[bucket_idx].resize(num_nodes);

    if (m_node_mass_fluid_y[bucket_idx].size() != num_nodes)
      m_node_mass_fluid_y[bucket_idx].resize(num_nodes);
    if (m_node_vel_fluid_y[bucket_idx].size() != num_nodes)
      m_node_vel_fluid_y[bucket_idx].resize(num_nodes);
    if (m_node_vol_fluid_y[bucket_idx].size() != num_nodes)
      m_node_vol_fluid_y[bucket_idx].resize(num_nodes);
    if (m_node_vol_pure_fluid_y[bucket_idx].size() != num_nodes)
      m_node_vol_pure_fluid_y[bucket_idx].resize(num_nodes);

    if (m_node_mass_fluid_z[bucket_idx].size() != num_nodes)
      m_node_mass_fluid_z[bucket_idx].resize(num_nodes);
    if (m_node_vel_fluid_z[bucket_idx].size() != num_nodes)
      m_node_vel_fluid_z[bucket_idx].resize(num_nodes);
    if (m_node_vol_fluid_z[bucket_idx].size() != num_nodes)
      m_node_vol_fluid_z[bucket_idx].resize(num_nodes);
    if (m_node_vol_pure_fluid_z[bucket_idx].size() != num_nodes)
      m_node_vol_pure_fluid_z[bucket_idx].resize(num_nodes);

    if (m_node_solid_phi[bucket_idx].size() != num_nodes)
      m_node_solid_phi[bucket_idx].resize(num_nodes);

    if (m_node_solid_vel_x[bucket_idx].size() != num_nodes)
      m_node_solid_vel_x[bucket_idx].resize(num_nodes);
    if (m_node_solid_vel_y[bucket_idx].size() != num_nodes)
      m_node_solid_vel_y[bucket_idx].resize(num_nodes);
    if (m_node_solid_vel_z[bucket_idx].size() != num_nodes)
      m_node_solid_vel_z[bucket_idx].resize(num_nodes);

    if (m_node_liquid_valid_x[bucket_idx].size() != num_nodes)
      m_node_liquid_valid_x[bucket_idx].resize(num_nodes);
    if (m_node_liquid_valid_y[bucket_idx].size() != num_nodes)
      m_node_liquid_valid_y[bucket_idx].resize(num_nodes);
    if (m_node_liquid_valid_z[bucket_idx].size() != num_nodes)
      m_node_liquid_valid_z[bucket_idx].resize(num_nodes);
  });
}

const std::vector<VectorXs>& TwoDScene::getNodeOrientationX() const {
  return m_node_orientation_x;
}

const std::vector<VectorXs>& TwoDScene::getNodeOrientationY() const {
  return m_node_orientation_y;
}

const std::vector<VectorXs>& TwoDScene::getNodeOrientationZ() const {
  return m_node_orientation_z;
}

const std::vector<VectorXs>& TwoDScene::getNodeShapeFactorX() const {
  return m_node_shape_factor_x;
}

const std::vector<VectorXs>& TwoDScene::getNodeShapeFactorY() const {
  return m_node_shape_factor_y;
}

const std::vector<VectorXs>& TwoDScene::getNodeShapeFactorZ() const {
  return m_node_shape_factor_z;
}

void TwoDScene::expandFluidNodesMarked(int layers) {
  auto check_bucket = [&](const Vector3i& bucket_handle,
                          const std::vector<unsigned char>& activated) {
    for (int t = -1; t <= 1; ++t)
      for (int s = -1; s <= 1; ++s)
        for (int r = -1; r <= 1; ++r) {
          if (t == 0 && s == 0 && r == 0) continue;

          Vector3i cur_bucket_handle = bucket_handle + Vector3i(r, s, t);
          if (cur_bucket_handle(0) < 0 ||
              cur_bucket_handle(0) >= m_particle_buckets.ni ||
              cur_bucket_handle(1) < 0 ||
              cur_bucket_handle(1) >= m_particle_buckets.nj ||
              cur_bucket_handle(2) < 0 ||
              cur_bucket_handle(2) >= m_particle_buckets.nk) {
            continue;
          }

          const int nbidx = m_particle_buckets.bucket_index(cur_bucket_handle);
          if (activated[nbidx]) {
            return true;
          }
        }

    return false;
  };

  for (int iLayer = 0; iLayer < layers; ++iLayer) {
    std::vector<unsigned char> activated_backup = m_bucket_activated;

    m_particle_buckets.for_each_bucket([&](int bucket_idx) {
      const Vector3i bucket_handle =
          m_particle_buckets.bucket_handle(bucket_idx);

      // ignore bucket already activated or no neighbor bucket activated
      if (activated_backup[bucket_idx] ||
          !check_bucket(bucket_handle, activated_backup))
        return;

      // Activate bucket that has activated neighbor buckets.
      m_bucket_activated[bucket_idx] = true;
    });
  }
}

/*!
 * resample nodes in the scene
 */
void TwoDScene::resampleNodes() {
  preAllocateNodes();

  auto particle_node_criteria = [this](int pidx) -> bool {
    return isSoft(pidx);
  };

  findNodes(m_particle_buckets, m_x, m_particle_nodes_x,
            Vector3s(0.0, 0.5, 0.5), particle_node_criteria);
  findNodes(m_particle_buckets, m_x, m_particle_nodes_y,
            Vector3s(0.5, 0.0, 0.5), particle_node_criteria);
  findNodes(m_particle_buckets, m_x, m_particle_nodes_z,
            Vector3s(0.5, 0.5, 0.0), particle_node_criteria);
  findNodes(m_particle_buckets, m_x, m_particle_nodes_solid_phi,
            Vector3s(0.0, 0.0, 0.0), particle_node_criteria);
  findNodes(m_particle_buckets, m_x, m_particle_nodes_p,
            Vector3s(0.5, 0.5, 0.5), particle_node_criteria);

  auto gauss_node_criteria = [this](int pidx) -> bool {
    return pidx < getNumEdges() + getNumFaces();
  };

  findNodes(m_gauss_buckets, m_x_gauss, m_gauss_nodes_x,
            Vector3s(0.0, 0.5, 0.5), gauss_node_criteria);
  findNodes(m_gauss_buckets, m_x_gauss, m_gauss_nodes_y,
            Vector3s(0.5, 0.0, 0.5), gauss_node_criteria);
  findNodes(m_gauss_buckets, m_x_gauss, m_gauss_nodes_z,
            Vector3s(0.5, 0.5, 0.0), gauss_node_criteria);

  expandFluidNodesMarked(1);

  // generate nodes in all activated buckets
  generateNodes();

  connectSolidPhiNodes();
  // connect pressure node to MAC nodes
  connectPressureNodes();

  markInsideOut();

  postAllocateNodes();
}

void TwoDScene::updateGaussManifoldSystem() {
  const int num_edges = m_edges.rows();

  threadutils::for_each(0, num_edges, [&](int i) {
    const auto& e = m_edges.row(i);
    m_fluid_vol_gauss(i) = (m_fluid_vol(e(0)) + m_fluid_vol(e(1))) * 0.5;
    m_v_gauss.segment<4>(i * 4) =
        (m_v.segment<4>(e(0) * 4) + m_v.segment<4>(e(1) * 4)) * 0.5;
    m_fluid_m_gauss.segment<4>(i * 4) =
        (m_fluid_m.segment<4>(e(0) * 4) + m_fluid_m.segment<4>(e(1) * 4)) * 0.5;
  });

  const int num_faces = m_faces.rows();
  threadutils::for_each(0, num_faces, [&](int i) {
    const auto& f = m_faces.row(i);

    const Vector3s& angle_frac = m_face_weights[i];

    m_fluid_vol_gauss(i + num_edges) = m_fluid_vol(f[0]) * angle_frac[0] +
                                       m_fluid_vol(f[1]) * angle_frac[1] +
                                       m_fluid_vol(f[2]) * angle_frac[2];
    m_v_gauss.segment<4>((i + num_edges) * 4) =
        m_v.segment<4>(f[0] * 4) * angle_frac[0] +
        m_v.segment<4>(f[1] * 4) * angle_frac[1] +
        m_v.segment<4>(f[2] * 4) * angle_frac[2];
    m_fluid_m_gauss.segment<4>((i + num_edges) * 4) =
        m_fluid_m.segment<4>(f[0] * 4) * angle_frac[0] +
        m_fluid_m.segment<4>(f[1] * 4) * angle_frac[1] +
        m_fluid_m.segment<4>(f[2] * 4) * angle_frac[2];
  });
}

void TwoDScene::updateGaussSystem(scalar dt) {
  const int num_edges = m_edges.rows();

  threadutils::for_each(0, num_edges, [&](int i) {
    const auto& e = m_edges.row(i);
    m_x_gauss.segment<4>(i * 4) =
        (m_x.segment<4>(e(0) * 4) + m_x.segment<4>(e(1) * 4)) * 0.5;
    m_v_gauss.segment<4>(i * 4) =
        (m_v.segment<4>(e(0) * 4) + m_v.segment<4>(e(1) * 4)) * 0.5;
    m_fluid_v_gauss.segment<4>(i * 4) =
        (m_fluid_v.segment<4>(e(0) * 4) + m_fluid_v.segment<4>(e(1) * 4)) * 0.5;
    m_fluid_vol_gauss(i) = (m_fluid_vol(e(0)) + m_fluid_vol(e(1))) * 0.5;
    m_fluid_m_gauss.segment<4>(i * 4) =
        (m_fluid_m.segment<4>(e(0) * 4) + m_fluid_m.segment<4>(e(1) * 4)) * 0.5;
  });

  // removed
  updateDeformationGradient(dt);
}

int TwoDScene::getNumElastoParticles() const {
  if (!m_fluids.size())
    return getNumParticles();
  else
    return m_fluids[0];
}

int TwoDScene::getNumSoftElastoParticles() const {
  return getNumElastoParticles() - getNumSurfels();
}

/*!
 * update plasticity for friction and sliding, see [Jiang et al. 2017]
 */
void TwoDScene::updatePlasticity(scalar dt) {
  const int num_edges = getNumEdges();
  // for curves
  threadutils::for_each(0, num_edges, [&](int pidx) {
    const Matrix3s& d_hat = m_d_gauss.block<3, 3>(pidx * 3, 0);
    Matrix3s Q, R;
    mathutils::QRDecompose<scalar, 3>(d_hat, Q, R);
    const scalar alpha = getFrictionAlpha(pidx);
    const scalar beta = getFrictionBeta(pidx);

    Eigen::JacobiSVD<Eigen::Matrix2d> svd(
        R.block<2, 2>(1, 1), Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector2d s = svd.singularValues();
    const Matrix2s U = svd.matrixU();
    const Matrix2s V = svd.matrixV();

    scalar ep1_hat = std::log(s(0));
    scalar ep2_hat = std::log(s(1));
    assert(ep1_hat >= ep2_hat);

    const scalar la = getLa(pidx) * getCollisionMultiplier(pidx);
    const scalar mu = getMu(pidx) * getCollisionMultiplier(pidx);

    Vector2s sigm_inv = Vector2s(1.0 / s(0), 1.0 / s(1));
    Vector2s lnsigm = Vector2s(ep1_hat, ep2_hat);

    if (ep1_hat + ep2_hat < 0) {
      Matrix2s ep = Matrix2s(lnsigm.asDiagonal());

      const scalar trep = ep.trace();
      Matrix2s eep = ep - trep * 0.5 * Matrix2s::Identity();

      scalar dgp = eep.norm() + (la + mu) / mu * trep * alpha;

      if (eep.norm() < 1e-20) {
        ep1_hat = 0.0;
        ep2_hat = 0.0;
      } else if (dgp > 0) {
        Matrix2s Hp = ep - dgp * eep / eep.norm();
        ep1_hat = Hp(0, 0);
        ep2_hat = Hp(1, 1);
      }
    } else {
      ep1_hat = 0.0;
      ep2_hat = 0.0;
    }

    s(0) = std::exp(ep1_hat);
    s(1) = std::exp(ep2_hat);

    sigm_inv = Vector2s(1.0 / s(0), 1.0 / s(1));
    lnsigm = Vector2s(ep1_hat, ep2_hat);

    R.block<2, 2>(1, 1) = U * Eigen::Matrix2d(s.asDiagonal()) * V.transpose();

    const scalar ff = mu * sqrt(R(0, 1) * R(0, 1) + R(0, 2) * R(0, 2));

    Matrix2s tmp =
        Matrix2s(sigm_inv.asDiagonal()) * Matrix2s(lnsigm.asDiagonal());
    const scalar fn =
        (2.0 * mu * tmp + la * lnsigm.sum() * Matrix2s(sigm_inv.asDiagonal()))
            .norm() *
        0.5;

    if (ff > 0.0 && ff > fn * beta) {
      R(0, 1) *= std::min(1.0, beta * fn / ff);
      R(0, 2) *= std::min(1.0, beta * fn / ff);
    }

    Matrix3s dhat = Q * R;

    m_Fe_gauss.block<3, 3>(3 * pidx, 0) =
        dhat * m_D_inv_gauss.block<3, 3>(3 * pidx, 0);
    m_d_gauss.block<3, 3>(3 * pidx, 0) = dhat;
  });

  // for cloth and surfels (removed)
}

Eigen::Quaternion<scalar>& TwoDScene::getGroupRotation(int group_idx) {
  return m_group_rot[group_idx];
}

Vector3s& TwoDScene::getGroupTranslation(int group_idx) {
  return m_group_pos[group_idx];
}

Eigen::Quaternion<scalar>& TwoDScene::getPrevGroupRotation(int group_idx) {
  return m_group_prev_rot[group_idx];
}

Vector3s& TwoDScene::getPrevGroupTranslation(int group_idx) {
  return m_group_prev_pos[group_idx];
}

void TwoDScene::resizeGroups(int num_group) {
  m_group_rot.resize(num_group);
  m_group_prev_rot.resize(num_group);
  m_group_pos.resize(num_group);
  m_group_prev_pos.resize(num_group);

  threadutils::for_each(0, num_group, [&](int i) {
    m_group_rot[i].setIdentity();
    m_group_prev_rot[i].setIdentity();
    m_group_pos[i].setZero();
    m_group_prev_pos[i].setZero();
  });

  m_shooting_vol_accum.resize(num_group, 0.0);
}

void TwoDScene::initGroupPos() {
  const int num_group = m_group_pos.size();

  for (int i = 0; i < num_group; ++i) {
    Vector3s center = Vector3s::Zero();
    m_group_prev_pos[i] = m_group_pos[i] = center;
  }
}

SimInfo& TwoDScene::getSimInfo() { return m_sim_info; }

const SimInfo& TwoDScene::getSimInfo() const { return m_sim_info; }

const std::vector<int>& TwoDScene::getFluidIndices() const { return m_fluids; }

std::vector<int>& TwoDScene::getFluidIndices() { return m_fluids; }

int TwoDScene::getNumFluidParticles() const { return (int)m_fluids.size(); }

scalar TwoDScene::computePhi(
    const Vector3s& pos)
    const {
  scalar min_phi = 3.0 * m_bucket_size;
  return min_phi;
}

scalar TwoDScene::computePhiVel(
    const Vector3s& pos, Vector3s& vel)
    const {
  scalar min_phi = 3.0 * m_bucket_size;
  Vector3s min_vel = Vector3s::Zero();
  vel = min_vel;

  return min_phi;
}

const std::vector<int>& TwoDScene::getParticleToSurfels() const {
  return m_particle_to_surfel;
}

const MatrixXs& TwoDScene::getGaussNormal() const { return m_norm_gauss; }

MatrixXs& TwoDScene::getGaussNormal() { return m_norm_gauss; }

/*!
 * update deformation gradient stored on face/edge
 */
void TwoDScene::updateDeformationGradient(scalar dt) {
  // updating deformation gradient
  const int num_edges = m_edges.rows();
  const scalar invD = getInverseDCoeff();

  threadutils::for_each(0, num_edges, [&](int pidx) {
    const auto& e = m_edges.row(pidx);
    auto& indices_x = m_gauss_nodes_x[pidx];
    auto& indices_y = m_gauss_nodes_y[pidx];
    auto& indices_z = m_gauss_nodes_z[pidx];

    auto& weights = m_gauss_weights[pidx];

    const Vector3s& pos = m_x_gauss.segment<3>(pidx * 4);

    Matrix3s gradx_hat;
    gradx_hat.setZero();

    for (int i = 0; i < indices_x.rows(); i++) {
      const int node_bucket_idx = indices_x(i, 0);
      const int node_idx = indices_x(i, 1);
      if (!m_bucket_activated[node_bucket_idx])
        continue;  // this shouldn't be touched

      const scalar& nv = m_node_vel_x[node_bucket_idx](node_idx);
      Vector3s np = getNodePosX(node_bucket_idx, node_idx);
      gradx_hat.block<1, 3>(0, 0) +=
          nv * weights(i, 0) * (np - pos).transpose() * invD;
      //      std::cout << "xi: " << i << ", " << g_grad_x.row(i) << std::endl;
    }

    for (int i = 0; i < indices_y.rows(); i++) {
      const int node_bucket_idx = indices_y(i, 0);
      const int node_idx = indices_y(i, 1);
      if (!m_bucket_activated[node_bucket_idx]) continue;

      const scalar& nv = m_node_vel_y[node_bucket_idx](node_idx);
      Vector3s np = getNodePosY(node_bucket_idx, node_idx);
      gradx_hat.block<1, 3>(1, 0) +=
          nv * weights(i, 1) * (np - pos).transpose() * invD;
      //      std::cout << "yi: " << i << ", " << g_grad_y.row(i) << std::endl;
    }

    for (int i = 0; i < indices_z.rows(); i++) {
      const int node_bucket_idx = indices_z(i, 0);
      const int node_idx = indices_z(i, 1);
      if (!m_bucket_activated[node_bucket_idx]) continue;

      const scalar& nv = m_node_vel_z[node_bucket_idx](node_idx);
      Vector3s np = getNodePosZ(node_bucket_idx, node_idx);
      gradx_hat.block<1, 3>(2, 0) +=
          nv * weights(i, 2) * (np - pos).transpose() * invD;
      //      std::cout << "zi: " << i << ", " << g_grad_z.row(i) << std::endl;
    }

    //        std::cout<<"gradx_hat: \n"<<gradx_hat<<std::endl;
    Matrix3s d_hat;

    Vector3s tangent = (m_x.segment<3>(e(1) * 4) - m_x.segment<3>(e(0) * 4));

    d_hat.block<3, 1>(0, 0) = tangent;
    d_hat.block<3, 2>(0, 1) = (Matrix3s::Identity() + gradx_hat * dt +
                               0.5 * gradx_hat * gradx_hat * (dt * dt)) *
                              m_d_gauss.block<3, 2>(pidx * 3, 1);

    m_d_gauss.block<3, 3>(pidx * 3, 0) = d_hat;
    m_Fe_gauss.block<3, 3>(pidx * 3, 0) =
        d_hat * m_D_inv_gauss.block<3, 3>(pidx * 3, 0);

    // update volume & volume fraction (removed)

    Matrix3s Q, R;
    mathutils::QRDecompose<scalar, 3>(d_hat, Q, R);

    m_norm_gauss.block<3, 3>(pidx * 3, 0) = Q;

    //        std::cout<<"m_d_gauss: \n"<<m_d_gauss.block<3, 3>(pidx*3
    //        ,0)<<std::endl; std::cout<<"gradx_hat: \n" <<gradx_hat<<std::endl;
    //    std::cout << std::endl;
  });

  // removed
}

/*!
 * update solid volume fraction
 */
void TwoDScene::updateSolidVolumeFraction() {
  const int num_soft_elasto = getNumSoftElastoParticles();
  const int num_edges = getNumEdges();

  threadutils::for_each(0, num_soft_elasto, [&](int pidx) {
    const std::vector<int>& edges = m_particle_to_edge[pidx];
    const std::vector<std::pair<int, scalar> >& faces =
        m_particle_to_face[pidx];

    scalar J = 0.0;
    scalar w = 0.0;

    for (int eidx : edges) {
      J += m_vol_gauss(eidx) * 0.5;
      w += m_rest_vol_gauss(eidx) * 0.5;
    }

    for (auto& p : faces) {
      const int gidx = p.first + num_edges;
      J += m_vol_gauss(gidx) * p.second;
      w += m_rest_vol_gauss(gidx) * p.second;
    }

    if (J > 1e-20 && w > 1e-20) {
      m_volume_fraction[pidx] = m_rest_volume_fraction[pidx] / J * w;
      m_vol[pidx] = m_rest_vol[pidx] * J / w;
    }
  });

  assert(!std::isnan(m_volume_fraction.sum()));
}

scalar TwoDScene::getCellSize() const {
  return m_bucket_size / (scalar)m_num_nodes;
}

scalar TwoDScene::getInverseDCoeff() const {
  return mathutils::inverse_D_coeff(getCellSize(), m_kernel_order);
}

void TwoDScene::setVolumeFraction(int particle, const scalar& vol_frac) {
  m_volume_fraction[particle] = m_rest_volume_fraction[particle] = vol_frac;
}

/*!
 * do DDA analysis for particle distribution
 */
void TwoDScene::computeDDA() {
  const int num_elasto = getNumElastoParticles();

  const int num_bin = 256;

  const scalar max_dist = getBucketLength();

  const scalar inc_dist = max_dist / (scalar)num_bin;

  const int num_buckets = getNumBuckets();

  std::vector<std::vector<int> > bucket_bins(num_buckets,
                                             std::vector<int>(num_bin, 0));

  m_particle_buckets.for_each_bucket_particles([&](int pidx, int bucket_idx) {
    if (pidx < num_elasto) return;

    std::vector<int>& bbin = bucket_bins[bucket_idx];

    m_particle_buckets.loop_neighbor_bucket_particles(bucket_idx, [&](int npidx,
                                                                      int) {
      if (npidx == pidx || npidx < num_elasto) return false;

      const scalar dist =
          (m_x.segment<3>(npidx * 4) - m_x.segment<3>(pidx * 4)).norm();
      const int bin_idx = (int)((dist - 0.5 * inc_dist) / max_dist * num_bin);
      if (bin_idx < 0 || bin_idx >= num_bin) return false;

      bbin[bin_idx] += 8;
      bbin[std::max(0, bin_idx - 1)] += 5;
      bbin[std::min(num_bin - 1, bin_idx + 1)] += 5;

      return false;
    });
  });

  std::vector<scalar> final_bins(num_bin, 0.0);

  scalar sum_bins = 0.0;
  for (int j = 0; j < num_bin; ++j) {
    for (int i = 0; i < num_buckets; ++i) {
      final_bins[j] += (scalar)bucket_bins[i][j];
    }
    const scalar dist = (scalar)j / (scalar)num_bin * max_dist + inc_dist * 0.5;
    if (dist > 0.0) final_bins[j] /= (dist * dist);

    sum_bins += final_bins[j];
  }

  if (sum_bins == 0.0) sum_bins = 1.0;

  for (int j = 0; j < num_bin; ++j) {
    final_bins[j] /= sum_bins;
  }

  std::cout << "[DDA Analysis]" << std::endl;
  std::cout << "-----------------------------------" << std::endl;
  for (int j = 0; j < num_bin; ++j) {
    const scalar dist = (scalar)j / (scalar)num_bin * max_dist + inc_dist * 0.5;
    std::cout << dist << ", " << final_bins[j] << std::endl;
  }
  std::cout << "-----------------------------------" << std::endl;
}

/*!
 * use MLS-MPM to map particle/vertices onto nodes, see [Hu et al. 2018]
 */
void TwoDScene::mapParticleNodesAPIC() {
  const scalar dx = getCellSize();
  const scalar dV = dx * dx * dx;
  //    std::cout << "FVb: " << m_fluid_v << std::endl;
  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    m_node_mass_x[bucket_idx].setZero();
    m_node_vel_x[bucket_idx].setZero();
    m_node_vol_x[bucket_idx].setZero();

    m_node_mass_y[bucket_idx].setZero();
    m_node_vel_y[bucket_idx].setZero();
    m_node_vol_y[bucket_idx].setZero();

    m_node_mass_z[bucket_idx].setZero();
    m_node_vel_z[bucket_idx].setZero();
    m_node_vol_z[bucket_idx].setZero();

    m_node_mass_fluid_x[bucket_idx].setZero();
    m_node_vel_fluid_x[bucket_idx].setZero();
    m_node_vol_fluid_x[bucket_idx].setZero();

    m_node_mass_fluid_y[bucket_idx].setZero();
    m_node_vel_fluid_y[bucket_idx].setZero();
    m_node_vol_fluid_y[bucket_idx].setZero();

    m_node_mass_fluid_z[bucket_idx].setZero();
    m_node_vel_fluid_z[bucket_idx].setZero();
    m_node_vol_fluid_z[bucket_idx].setZero();

    m_node_psi_x[bucket_idx].setZero();
    m_node_sat_x[bucket_idx].setZero();

    m_node_psi_y[bucket_idx].setZero();
    m_node_sat_y[bucket_idx].setZero();

    m_node_psi_z[bucket_idx].setZero();
    m_node_sat_z[bucket_idx].setZero();

    const auto& bucket_node_particles_x = m_node_particles_x[bucket_idx];
    const auto& bucket_node_particles_y = m_node_particles_y[bucket_idx];
    const auto& bucket_node_particles_z = m_node_particles_z[bucket_idx];

    const int num_nodes = getNumNodes(bucket_idx);

    for (int i = 0; i < num_nodes; ++i) {
      const auto& node_particles_x = bucket_node_particles_x[i];
      const Vector3s& np = getNodePosX(bucket_idx, i);

      scalar p = 0.0;
      scalar mass = 0.0;
      scalar vol_solid = 0.0;

      scalar p_fluid = 0.0;
      scalar mass_fluid = 0.0;
      scalar vol_fluid = 0.0;

      scalar vol_fluid_elasto = 0.0;
      scalar shape_factor = 0.0;
      scalar shape_factor_rw = 0.0;
      Vector3s orientation = Vector3s::Zero();

      for (auto& pair : node_particles_x) {
        const int pidx = pair.first;

        auto& weights = m_particle_weights[pidx];
        const scalar& vol = m_rest_vol(pidx);
        const Vector3s& m = m_m.segment<3>(pidx * 4);
        const scalar& fvol = m_fluid_vol(pidx);
        const Vector3s& fm = m_fluid_m.segment<3>(pidx * 4);
        const Vector3s& v = m_v.segment<3>(pidx * 4);
        const Vector3s& fluidv = m_fluid_v.segment<3>(pidx * 4);
        const Vector3s& pos = m_x.segment<3>(pidx * 4);
        const Matrix3s& B = m_B.block<3, 3>(pidx * 3, 0);
        const Matrix3s& fB = m_fB.block<3, 3>(pidx * 3, 0);

        if (!isFluid(pidx)) {
          const scalar vel = v(0) + B.row(0).dot(np - pos);
          p += vel * (m(0) + fm(0)) * weights(pair.second, 0);
          mass += (m(0) + fm(0)) * weights(pair.second, 0);

          if (m_particle_to_surfel[pidx] < 0) {
            vol_solid +=
                vol * m_rest_volume_fraction(pidx) * weights(pair.second, 0);
            vol_fluid_elasto += fvol * weights(pair.second, 0);
            shape_factor += m_shape_factor(pidx) * weights(pair.second, 0);
            shape_factor_rw += weights(pair.second, 0);
            orientation +=
                m_orientation.segment<3>(pidx * 3) * weights(pair.second, 0);
          }
        } else {
          const scalar vel = fluidv(0) + fB.row(0).dot(np - pos);
          p_fluid += vel * fm(0) * weights(pair.second, 0);
          mass_fluid += fm(0) * weights(pair.second, 0);
          vol_fluid += fvol * weights(pair.second, 0);
        }
      }

      if (mass > 1e-20) {
        m_node_vel_x[bucket_idx](i) = p / mass;
      }

      if (mass_fluid > 1e-20) {
        m_node_vel_fluid_x[bucket_idx](i) = p_fluid / mass_fluid;
      }

      if (shape_factor_rw > 1e-20) {
        shape_factor /= shape_factor_rw;
      }

      m_node_mass_x[bucket_idx](i) = mass;
      m_node_vol_x[bucket_idx](i) = vol_solid + vol_fluid_elasto;

      m_node_mass_fluid_x[bucket_idx](i) = mass_fluid;
      m_node_vol_fluid_x[bucket_idx](i) = vol_fluid;

      m_node_psi_x[bucket_idx](i) = mathutils::clamp(vol_solid / dV, 0.0, 1.0);
      m_node_sat_x[bucket_idx](i) = mathutils::clamp(
          (vol_fluid + vol_fluid_elasto) / std::max(1e-20, dV - vol_solid), 0.0,
          1.0);

      const scalar lo = orientation.norm();
      if (lo > 1e-20) {
        orientation /= lo;
      }
      m_node_orientation_x[bucket_idx].segment<3>(i * 3) = orientation;
      m_node_shape_factor_x[bucket_idx](i) = shape_factor;
    }

    assert(!std::isnan(m_node_vel_x[bucket_idx].sum()));
    assert(!std::isnan(m_node_mass_fluid_x[bucket_idx].sum()));
    assert(!std::isnan(m_node_vol_fluid_x[bucket_idx].sum()));

    for (int i = 0; i < num_nodes; ++i) {
      const auto& node_particles_y = bucket_node_particles_y[i];
      const Vector3s& np = getNodePosY(bucket_idx, i);

      scalar p = 0.0;
      scalar mass = 0.0;
      scalar vol_solid = 0.0;

      scalar p_fluid = 0.0;
      scalar mass_fluid = 0.0;
      scalar vol_fluid = 0.0;

      scalar vol_fluid_elasto = 0.0;
      scalar shape_factor = 0.0;
      scalar shape_factor_rw = 0.0;
      Vector3s orientation = Vector3s::Zero();

      for (auto& pair : node_particles_y) {
        const int pidx = pair.first;

        auto& weights = m_particle_weights[pidx];
        const scalar& vol = m_rest_vol(pidx);
        const Vector3s& m = m_m.segment<3>(pidx * 4);
        const scalar& fvol = m_fluid_vol(pidx);
        const Vector3s& fm = m_fluid_m.segment<3>(pidx * 4);
        const Vector3s& v = m_v.segment<3>(pidx * 4);
        const Vector3s& fluidv = m_fluid_v.segment<3>(pidx * 4);
        const Vector3s& pos = m_x.segment<3>(pidx * 4);
        const Matrix3s& B = m_B.block<3, 3>(pidx * 3, 0);
        const Matrix3s& fB = m_fB.block<3, 3>(pidx * 3, 0);

        if (!isFluid(pidx)) {
          const scalar vel = v(1) + B.row(1).dot(np - pos);
          p += vel * (m(1) + fm(1)) * weights(pair.second, 1);
          mass += (m(1) + fm(1)) * weights(pair.second, 1);

          if (m_particle_to_surfel[pidx] < 0) {
            vol_solid +=
                vol * m_rest_volume_fraction(pidx) * weights(pair.second, 1);
            vol_fluid_elasto += fvol * weights(pair.second, 1);
            shape_factor += m_shape_factor(pidx) * weights(pair.second, 1);
            shape_factor_rw += weights(pair.second, 1);
            orientation +=
                m_orientation.segment<3>(pidx * 3) * weights(pair.second, 1);
          }
        } else {
          const scalar vel = fluidv(1) + fB.row(1).dot(np - pos);
          p_fluid += vel * fm(1) * weights(pair.second, 1);
          mass_fluid += fm(1) * weights(pair.second, 1);
          vol_fluid += fvol * weights(pair.second, 1);
        }
      }

      if (mass > 1e-20) {
        m_node_vel_y[bucket_idx](i) = p / mass;
      }

      if (mass_fluid > 1e-20) {
        m_node_vel_fluid_y[bucket_idx](i) = p_fluid / mass_fluid;
      }

      if (shape_factor_rw > 1e-20) {
        shape_factor /= shape_factor_rw;
      }

      m_node_mass_y[bucket_idx](i) = mass;
      m_node_vol_y[bucket_idx](i) = vol_solid + vol_fluid_elasto;

      m_node_mass_fluid_y[bucket_idx](i) = mass_fluid;
      m_node_vol_fluid_y[bucket_idx](i) = vol_fluid;

      m_node_psi_y[bucket_idx](i) = mathutils::clamp(vol_solid / dV, 0.0, 1.0);
      m_node_sat_y[bucket_idx](i) = mathutils::clamp(
          (vol_fluid + vol_fluid_elasto) / std::max(1e-20, dV - vol_solid), 0.0,
          1.0);

      const scalar lo = orientation.norm();
      if (lo > 1e-20) {
        orientation /= lo;
      }
      m_node_orientation_y[bucket_idx].segment<3>(i * 3) = orientation;
      m_node_shape_factor_y[bucket_idx](i) = shape_factor;
    }

    assert(!std::isnan(m_node_vel_y[bucket_idx].sum()));
    assert(!std::isnan(m_node_mass_fluid_y[bucket_idx].sum()));
    assert(!std::isnan(m_node_vol_fluid_y[bucket_idx].sum()));

    for (int i = 0; i < num_nodes; ++i) {
      const auto& node_particles_z = bucket_node_particles_z[i];
      const Vector3s& np = getNodePosZ(bucket_idx, i);

      scalar p = 0.0;
      scalar mass = 0.0;
      scalar vol_solid = 0.0;

      scalar p_fluid = 0.0;
      scalar mass_fluid = 0.0;
      scalar vol_fluid = 0.0;

      scalar vol_fluid_elasto = 0.0;
      scalar shape_factor = 0.0;
      scalar shape_factor_rw = 0.0;
      Vector3s orientation = Vector3s::Zero();

      for (auto& pair : node_particles_z) {
        const int pidx = pair.first;

        auto& weights = m_particle_weights[pidx];
        const scalar& vol = m_rest_vol(pidx);
        const Vector3s& m = m_m.segment<3>(pidx * 4);
        const scalar& fvol = m_fluid_vol(pidx);
        const Vector3s& fm = m_fluid_m.segment<3>(pidx * 4);
        const Vector3s& v = m_v.segment<3>(pidx * 4);
        const Vector3s& fluidv = m_fluid_v.segment<3>(pidx * 4);
        const Vector3s& pos = m_x.segment<3>(pidx * 4);
        const Matrix3s& B = m_B.block<3, 3>(pidx * 3, 0);
        const Matrix3s& fB = m_fB.block<3, 3>(pidx * 3, 0);

        if (!isFluid(pidx)) {
          const scalar vel = v(2) + B.row(2).dot(np - pos);
          p += vel * (m(2) + fm(2)) * weights(pair.second, 2);
          mass += (m(2) + fm(2)) * weights(pair.second, 2);

          if (m_particle_to_surfel[pidx] < 0) {
            vol_solid +=
                vol * m_rest_volume_fraction(pidx) * weights(pair.second, 2);
            vol_fluid_elasto += fvol * weights(pair.second, 2);
            shape_factor += m_shape_factor(pidx) * weights(pair.second, 2);
            shape_factor_rw += weights(pair.second, 2);
            orientation +=
                m_orientation.segment<3>(pidx * 3) * weights(pair.second, 2);
          }
        } else {
          const scalar vel = fluidv(2) + fB.row(2).dot(np - pos);
          p_fluid += vel * fm(2) * weights(pair.second, 2);
          mass_fluid += fm(2) * weights(pair.second, 2);
          vol_fluid += fvol * weights(pair.second, 2);
        }
      }

      if (mass > 1e-20) {
        m_node_vel_z[bucket_idx](i) = p / mass;
      }

      if (mass_fluid > 1e-20) {
        m_node_vel_fluid_z[bucket_idx](i) = p_fluid / mass_fluid;
      }

      if (shape_factor_rw > 1e-20) {
        shape_factor /= shape_factor_rw;
      }

      m_node_mass_z[bucket_idx](i) = mass;
      m_node_vol_z[bucket_idx](i) = vol_solid + vol_fluid_elasto;

      m_node_mass_fluid_z[bucket_idx](i) = mass_fluid;
      m_node_vol_fluid_z[bucket_idx](i) = vol_fluid;

      m_node_psi_z[bucket_idx](i) = mathutils::clamp(vol_solid / dV, 0.0, 1.0);
      m_node_sat_z[bucket_idx](i) = mathutils::clamp(
          (vol_fluid + vol_fluid_elasto) / std::max(1e-20, dV - vol_solid), 0.0,
          1.0);

      const scalar lo = orientation.norm();
      if (lo > 1e-20) {
        orientation /= lo;
      }
      m_node_orientation_z[bucket_idx].segment<3>(i * 3) = orientation;
      m_node_shape_factor_z[bucket_idx](i) = shape_factor;
    }

    assert(!std::isnan(m_node_vel_z[bucket_idx].sum()));
    assert(!std::isnan(m_node_mass_fluid_z[bucket_idx].sum()));
    assert(!std::isnan(m_node_vol_fluid_z[bucket_idx].sum()));
  });
}

bool TwoDScene::isFluid(int pidx) const {
  return pidx >= getNumElastoParticles();
}

void TwoDScene::saveFluidVelocity() {
  m_node_vel_saved_fluid_x = m_node_vel_fluid_x;
  m_node_vel_saved_fluid_y = m_node_vel_fluid_y;
  m_node_vel_saved_fluid_z = m_node_vel_fluid_z;
}

void TwoDScene::saveParticleVelocity() { m_saved_v = m_v; }

/*!
 * map node variables back to particle and vertices, see [Jiang et al. 2015]
 */
void TwoDScene::mapNodeParticlesAPIC() {
  const int num_part = getNumParticles();

  const scalar invD = getInverseDCoeff();

  threadutils::for_each(0, num_part, [&](int pidx) {
    if (m_particle_to_surfel[pidx] >= 0 || isOutsideFluid(pidx)) return;

    auto& indices_x = m_particle_nodes_x[pidx];
    auto& indices_y = m_particle_nodes_y[pidx];
    auto& indices_z = m_particle_nodes_z[pidx];

    auto& weights = m_particle_weights[pidx];
    const Vector3s& pos = m_x.segment<3>(pidx * 4);

    m_v.segment<3>(pidx * 4).setZero();

    m_fluid_v.segment<4>(pidx * 4).setZero();

    m_B.block<3, 3>(pidx * 3, 0).setZero();
    m_fB.block<3, 3>(pidx * 3, 0).setZero();

    bool is_fluid = isFluid(pidx);

    for (int i = 0; i < indices_x.rows(); ++i) {
        const int node_bucket_idx = indices_x(i, 0);
        if (!m_bucket_activated[node_bucket_idx]) continue;

        const int node_idx = indices_x(i, 1);

        const Vector3s& np = getNodePosX(node_bucket_idx, node_idx);
        const scalar& nv = m_node_vel_x[node_bucket_idx](node_idx);

        m_v(pidx * 4 + 0) += nv * weights(i, 0);

        m_B.block<1, 3>(pidx * 3 + 0, 0) +=
            nv * weights(i, 0) * (np - pos).transpose() * invD;
    }

    assert(!std::isnan(m_v.segment<3>(pidx * 4).sum()));
    assert(!std::isnan(m_fluid_v.segment<3>(pidx * 4).sum()));

    for (int i = 0; i < indices_y.rows(); ++i) {
        const int node_bucket_idx = indices_y(i, 0);
        if (!m_bucket_activated[node_bucket_idx]) continue;

        const int node_idx = indices_y(i, 1);

        const Vector3s& np = getNodePosY(node_bucket_idx, node_idx);
        const scalar& nv = m_node_vel_y[node_bucket_idx](node_idx);

        m_v(pidx * 4 + 1) += nv * weights(i, 1);

        m_B.block<1, 3>(pidx * 3 + 1, 0) +=
            nv * weights(i, 1) * (np - pos).transpose() * invD;
    }

    assert(!std::isnan(m_v.segment<3>(pidx * 4).sum()));
    assert(!std::isnan(m_fluid_v.segment<3>(pidx * 4).sum()));

    for (int i = 0; i < indices_z.rows(); ++i) {
        const int node_bucket_idx = indices_z(i, 0);
        if (!m_bucket_activated[node_bucket_idx]) continue;

        const int node_idx = indices_z(i, 1);

        const Vector3s& np = getNodePosZ(node_bucket_idx, node_idx);
        const scalar& nv = m_node_vel_z[node_bucket_idx](node_idx);

        m_v(pidx * 4 + 2) += nv * weights(i, 2);

        m_B.block<1, 3>(pidx * 3 + 2, 0) +=
            nv * weights(i, 2) * (np - pos).transpose() * invD;
    }

    m_v.segment<4>(pidx * 4) *= m_sim_info.elasto_advect_coeff;

    m_B.block<3, 3>(pidx * 3, 0) =
        ((m_sim_info.elasto_flip_coeff +
            m_sim_info.elasto_flip_asym_coeff) *
            m_B.block<3, 3>(pidx * 3, 0) +
            (m_sim_info.elasto_flip_coeff -
                m_sim_info.elasto_flip_asym_coeff) *
            m_B.block<3, 3>(pidx * 3, 0).transpose()) *
        0.5;

    assert(!std::isnan(m_v.segment<3>(pidx * 4).sum()));
    assert(!std::isnan(m_fluid_v.segment<3>(pidx * 4).sum()));
  });
}

void TwoDScene::insertSolveGroup(const VectorXi& group) {
  m_solve_groups.push_back(group);
}

const std::vector<VectorXi>& TwoDScene::getSolveGroup() const {
  return m_solve_groups;
}

void TwoDScene::updateVelocityDifference() { m_dv = m_v - m_saved_v; }

Sorter& TwoDScene::getParticleBuckets() { return m_particle_buckets; }

const Sorter& TwoDScene::getParticleBuckets() const {
  return m_particle_buckets;
}

Sorter& TwoDScene::getGaussBuckets() { return m_gauss_buckets; }

const Sorter& TwoDScene::getGaussBuckets() const { return m_gauss_buckets; }

void TwoDScene::setPosition(int particle, const Vector3s& pos) {
  assert(particle >= 0);
  assert(particle < getNumParticles());
  //    m_x.segment<3>( getDof(particle) ) = pos;
  m_x.segment<3>(4 * particle) = pos;
}

void TwoDScene::setTheta(int particle, const scalar theta) {
  assert(particle >= 0);
  assert(particle < getNumParticles());
  m_x(4 * particle + 3) = theta;
}

void TwoDScene::setOmega(int particle, const scalar omega) {
  assert(particle >= 0);
  assert(particle < getNumParticles());
  m_v(4 * particle + 3) = omega;
}

void TwoDScene::setVelocity(int particle, const Vector3s& vel) {
  assert(particle >= 0);
  assert(particle < getNumParticles());

  m_v.segment<3>(4 * particle) = vel;
  m_B.block<3, 3>(3 * particle, 0).setZero();
}

void TwoDScene::setTipVerts(int particle, bool tipVerts) {
  m_is_strand_tip[particle] = tipVerts;
}

void TwoDScene::setVolume(int particle, const scalar& volume) {
  assert(particle >= 0);
  assert(particle < getNumParticles());
  m_rest_vol(particle) = m_vol(particle) = volume;
}

void TwoDScene::setFluidVolume(int particle, const scalar& volume) {
  assert(particle >= 0);
  assert(particle < getNumParticles());
  m_fluid_vol(particle) = volume;
}

void TwoDScene::setGroup(int particle, int group) {
  assert(particle >= 0);
  assert(particle < getNumParticles());
  m_particle_group[particle] = group;
}

void TwoDScene::setRadius(int particle, const scalar& radiusA,
                          const scalar& radiusB) {
  assert(particle >= 0);
  assert(particle < getNumParticles());
  m_radius(particle * 2 + 0) = radiusA;
  m_radius(particle * 2 + 1) = radiusB;
}

void TwoDScene::setMass(int particle, const scalar& mass,
                        const scalar& second_moments) {
  assert(particle >= 0);
  assert(particle < getNumParticles());

  m_m(4 * particle) = mass;
  m_m(4 * particle + 1) = mass;
  m_m(4 * particle + 2) = mass;
  m_m(4 * particle + 3) = second_moments;
}

void TwoDScene::setFluidMass(int particle, const scalar& mass,
                             const scalar& second_moments) {
  assert(particle >= 0);
  assert(particle < getNumParticles());

  m_fluid_m(4 * particle) = mass;
  m_fluid_m(4 * particle + 1) = mass;
  m_fluid_m(4 * particle + 2) = mass;
  m_fluid_m(4 * particle + 3) = second_moments;
}

void TwoDScene::updateShapeFactor() {
  const int num_elasto = getNumElastoParticles();
  threadutils::for_each(0, num_elasto, [&](int pidx) {
    if (m_particle_to_surfel[pidx] >= 0) {
      m_shape_factor(pidx) = 0.0;
    } else {
      const std::vector<int>& edges = m_particle_to_edge[pidx];
      const std::vector<std::pair<int, scalar> >& faces =
          m_particle_to_face[pidx];

      if (faces.size() == 0) {
        m_shape_factor(pidx) = 1.0;
      } else if (edges.size() == 0) {
        m_shape_factor(pidx) = 0.0;
      } else {
        scalar vol_edges = 0.0;

        for (int j : edges) {
          vol_edges += m_vol_gauss(j) * 0.5;
        }

        m_shape_factor(pidx) =
            mathutils::clamp(vol_edges / m_vol(pidx), 0.0, 1.0);
      }
    }
  });
}

void TwoDScene::updateOrientation() {
  const int num_elasto = getNumElastoParticles();
  const int num_edges = getNumEdges();
  threadutils::for_each(0, num_elasto, [&](int pidx) {
    if (m_particle_to_surfel[pidx] >= 0) {
      // removed
    } else {
      Vector3s ori = Vector3s::Zero();
      const std::vector<int>& edges = m_particle_to_edge[pidx];
      const std::vector<std::pair<int, scalar> >& faces =
          m_particle_to_face[pidx];

      for (int eidx : edges) {
        ori += (m_x.segment<3>(m_edges(eidx, 1) * 4) -
                m_x.segment<3>(m_edges(eidx, 0) * 4))
                   .normalized() *
               m_vol_gauss(eidx) * 0.5;
      }

      for (auto& p : faces) {
        // removed
      }

      m_orientation.segment<3>(pidx * 3) = ori.normalized();
    }
  });
}

void TwoDScene::setEdge(int idx, const std::pair<int, int>& edge) {
  m_edges(idx, 0) = edge.first;
  m_edges(idx, 1) = edge.second;
  m_edge_inv_mapping[idx] = Vector2i(m_particle_to_edge[edge.first].size(),
                                     m_particle_to_edge[edge.second].size());
  m_particle_to_edge[edge.first].push_back(idx);
  m_particle_to_edge[edge.second].push_back(idx);
}

void TwoDScene::setFace(int idx, const Vector3i& face) {
  m_faces.row(idx) = face.transpose();

  const Vector3s& x0 = m_rest_x.segment<3>(m_faces(idx, 0) * 4);
  const Vector3s& x1 = m_rest_x.segment<3>(m_faces(idx, 1) * 4);
  const Vector3s& x2 = m_rest_x.segment<3>(m_faces(idx, 2) * 4);

  Vector3s angle_frac = Vector3s::Constant(1.0 / 3.0);
  angle_frac[0] =
      atan2((x1 - x0).cross(x2 - x0).norm(), (x1 - x0).dot(x2 - x0)) / PI;
  angle_frac[1] =
      atan2((x0 - x1).cross(x2 - x1).norm(), (x0 - x1).dot(x2 - x1)) / PI;
  angle_frac[2] = 1. - angle_frac[0] - angle_frac[1];

  m_face_weights[idx] = angle_frac;

  m_face_inv_mapping[idx] = Vector3i(m_particle_to_face[face(0)].size(),
                                     m_particle_to_face[face(1)].size(),
                                     m_particle_to_face[face(2)].size());

  m_particle_to_face[face(0)].push_back(
      std::pair<int, scalar>(idx, angle_frac[0]));
  m_particle_to_face[face(1)].push_back(
      std::pair<int, scalar>(idx, angle_frac[1]));
  m_particle_to_face[face(2)].push_back(
      std::pair<int, scalar>(idx, angle_frac[2]));
}

void TwoDScene::setFixed(int particle, unsigned char fixed) {
  assert(particle >= 0);
  assert(particle < getNumParticles());

  m_fixed[particle] = fixed;
}

void TwoDScene::setTwist(int particle, bool twist) {
  assert(particle >= 0);
  assert(particle < getNumParticles());

  m_twist[particle] = twist;
}

unsigned char TwoDScene::isFixed(int particle) const {
  assert(particle >= 0);
  assert(particle < getNumParticles());

  return m_fixed[particle];
}

bool TwoDScene::isTwist(int particle) const {
  assert(particle >= 0);
  assert(particle < getNumParticles());

  return m_twist[particle];
}

VectorXs TwoDScene::getPosition(int particle) {
  assert(particle >= 0);
  // assert( particle < getNumParticles() );
  assert(getDof(particle) < m_x.size());

  return m_x.segment<3>(getDof(particle));
}

int TwoDScene::getDof(int particle) const { return particle * 4; }

void TwoDScene::insertElasticParameters(
    const std::shared_ptr<ElasticParameters>& newparams) {
  m_strandParameters.push_back(newparams);
}

std::shared_ptr<ElasticParameters>& TwoDScene::getElasticParameters(
    const int index) {
  assert(0 <= index);
  assert(index < m_strandParameters.size());
  return m_strandParameters[index];
}

void TwoDScene::setEdgeRestLength(int idx, const scalar& l0) {
  m_edge_rest_length(idx) = l0;

  m_particle_rest_length(m_edges(idx, 0)) += l0 * 0.5;
  m_particle_rest_length(m_edges(idx, 1)) += l0 * 0.5;

  m_particle_rest_area(m_edges(idx, 0)) +=
      l0 *
      mathutils::perimeter(m_radius(m_edges(idx, 0) * 2 + 0),
                           m_radius(m_edges(idx, 0) * 2 + 1)) *
      0.5;
  m_particle_rest_area(m_edges(idx, 1)) +=
      l0 *
      mathutils::perimeter(m_radius(m_edges(idx, 1) * 2 + 0),
                           m_radius(m_edges(idx, 1) * 2 + 1)) *
      0.5;
}

void TwoDScene::setFaceRestArea(int idx, const scalar& a0) {
  m_face_rest_area(idx) = a0;

  for (unsigned int n = 0; n < 3; n++)
    m_particle_rest_area(m_faces(idx, n)) += a0 / 3.0;
}

void TwoDScene::updateRestPos() { m_rest_x = m_x; }

const VectorXs& TwoDScene::getRestPos() const { return m_rest_x; }

VectorXs& TwoDScene::getRestPos() { return m_rest_x; }

Vector3s TwoDScene::getTwistDir(int particle) const {
  auto& edges = m_particle_to_edge[particle];

  Vector3s dir = Vector3s::Zero();

  for (int eidx : edges) {
    dir += m_x.segment<3>(m_edges(eidx, 1) * 4) -
           m_x.segment<3>(m_edges(eidx, 0) * 4);
  }

  return dir.normalized();
}

Vector3s TwoDScene::getRestTwistDir(int particle) const {
  auto& edges = m_particle_to_edge[particle];

  Vector3s dir = Vector3s::Zero();

  for (int eidx : edges) {
    dir += m_rest_x.segment<3>(m_edges(eidx, 1) * 4) -
           m_rest_x.segment<3>(m_edges(eidx, 0) * 4);
  }

  return dir.normalized();
}

scalar TwoDScene::getParticleRestArea(int idx) const {
  return m_particle_rest_area(idx);
}

scalar TwoDScene::getParticleRestLength(int idx) const {
  return m_particle_rest_length(idx);
}

void TwoDScene::clearEdges() { m_edges.resize(0, 3); }

const VectorXs& TwoDScene::getFaceRestArea() const { return m_face_rest_area; }

const VectorXs& TwoDScene::getEdgeRestLength() const {
  return m_edge_rest_length;
}

const MatrixXi& TwoDScene::getFaces() const { return m_faces; }

const MatrixXi& TwoDScene::getEdges() const { return m_edges; }

const std::vector<int>& TwoDScene::getSurfels() const { return m_surfels; }

const std::vector<std::shared_ptr<AttachForce> >& TwoDScene::getAttachForces()
    const {
  return m_attach_forces;
}

const Vector2iT TwoDScene::getEdge(int edg) const {
  assert(edg >= 0);
  assert(edg < (int)m_edges.rows());

  return m_edges.row(edg);
}

void TwoDScene::insertScript(const std::shared_ptr<Script>& script) {
  m_scripts.push_back(script);
}

void TwoDScene::insertForce(const std::shared_ptr<Force>& newforce) {
  m_forces.push_back(newforce);
}

void TwoDScene::insertStrandForce(const std::shared_ptr<StrandForce>& newforce) {
  m_strands.push_back(newforce);
}

void TwoDScene::setStrandForceScript(int strand_idx, int node, scalar t)
{
  m_strands[strand_idx]->AddScriptedForce(node, t);
}

scalar TwoDScene::computeKineticEnergy() const {
  scalar T = 0.0;
  for (int i = 0; i < getNumParticles(); ++i) {
    T += m_m(4 * i) * m_v.segment<3>(4 * i).dot(m_v.segment<3>(4 * i));
  }
  return 0.5 * T;
}

scalar TwoDScene::computePotentialEnergy() const {
  scalar U = 0.0;
  for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i)
    m_forces[i]->addEnergyToTotal(m_x, m_v, m_m, m_volume_fraction,
                                  m_sim_info.lambda, U);
  return U;
}

void TwoDScene::postcompute(VectorXs& v, const scalar& dt) {
  threadutils::for_each(0, (int)m_forces.size(),
                        [&](int f) { m_forces[f]->postCompute(v, dt); });
}

void TwoDScene::precompute() {
  threadutils::for_each(0, (int)m_forces.size(),
                        [&](int f) { m_forces[f]->preCompute(); });
}

void TwoDScene::updateStartState() {
  threadutils::for_each(0, (int)m_forces.size(),
                        [&](int f) { m_forces[f]->updateStartState(); });
}

scalar TwoDScene::computeTotalEnergy() const {
  return computeKineticEnergy() + computePotentialEnergy();
}

const VectorXs& TwoDScene::getVolumeFraction() const {
  return m_volume_fraction;
}

VectorXs& TwoDScene::getVolumeFraction() { return m_volume_fraction; }

/*!
 * update differential operators on rods and meshes
 */
void TwoDScene::updateManifoldOperators() {
  const int num_edges = m_edges.rows();
  const int num_triangles = m_faces.rows();
  const int num_surfels = m_surfels.size();

  threadutils::for_each(0, num_edges, [&](int i) {
    const auto& e = m_edges.row(i);

    m_grad_gauss.block<3, 3>(i * 3, 0).setZero();
    Vector3s ev = m_x.segment<3>(e(1) * 4) - m_x.segment<3>(e(0) * 4);
    scalar l2ev = ev.squaredNorm();
    if (l2ev > 1e-20) ev /= l2ev;
    m_grad_gauss.block<3, 1>(i * 3, 0) = -ev;
    m_grad_gauss.block<3, 1>(i * 3, 1) = ev;
  });

  threadutils::for_each(0, num_triangles, [&](int i) {
    const auto& f = m_faces.row(i);

    Matrix3s g;
    mathutils::grad_triangle(m_x.segment<3>(f[0] * 4), m_x.segment<3>(f[1] * 4),
                             m_x.segment<3>(f[2] * 4), g);
    m_grad_gauss.block<3, 3>((i + num_edges) * 3, 0) = g;
  });

  threadutils::for_each(0, num_surfels, [&](int i) {
    int gidx = i + num_edges + num_triangles;

    m_grad_gauss.block<3, 3>(gidx * 3, 0).setZero();
  });

  updateParticleDiv();
}

void TwoDScene::updateGaussAccel() {
  const int num_edges = m_edges.rows();
  const int num_triangles = m_faces.rows();
  const int num_surfels = m_surfels.size();

  // update gauss dv
  threadutils::for_each(0, num_edges, [&](int i) {
    const auto& e = m_edges.row(i);
    m_dv_gauss.segment<4>(i * 4) =
        (m_dv.segment<4>(e(0) * 4) + m_dv.segment<4>(e(1) * 4)) * 0.5;
  });

  // remove
}

const VectorXs& TwoDScene::getGaussDV() const { return m_dv_gauss; }

VectorXs& TwoDScene::getGaussDV() { return m_dv_gauss; }

const VectorXs& TwoDScene::getGaussFluidM() const { return m_fluid_m_gauss; }

VectorXs& TwoDScene::getGaussFluidM() { return m_fluid_m_gauss; }

const VectorXs& TwoDScene::getGaussFluidVol() const {
  return m_fluid_vol_gauss;
}

VectorXs& TwoDScene::getGaussFluidVol() { return m_fluid_vol_gauss; }

const VectorXs& TwoDScene::getGaussVolumeFraction() const {
  return m_volume_fraction_gauss;
}

VectorXs& TwoDScene::getGaussVolumeFraction() {
  return m_volume_fraction_gauss;
}

const std::vector<VectorXs>& TwoDScene::getParticleDiv() const { return m_div; }

const std::vector<int>& TwoDScene::getParticleEdges(int pidx) const {
  return m_particle_to_edge[pidx];
}

const std::vector<std::pair<int, scalar> >& TwoDScene::getParticleFaces(
    int pidx) const {
  return m_particle_to_face[pidx];
}

/*!
 * calculate gradient of velocity on manifold
 */
void TwoDScene::accumulateManifoldFluidGradU(VectorXs& F) {
  const int ndof = getNumParticles() * 4;

  VectorXs F_full(ndof);
  F_full.setZero();

  for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i) {
    if (m_forces[i]->flag() & 2)
      m_forces[i]->addGradEToTotal(m_x, m_fluid_v, m_fluid_m, m_volume_fraction,
                                   m_sim_info.lambda, F_full);
  }

  F_full *= -1.0;

  // project to manifold
  const int num_edges = getNumEdges();
  const int num_faces = getNumFaces();

  threadutils::for_each(0, num_edges, [&](int eidx) {
    const auto& e = m_edges.row(eidx);
    Vector3s fu =
        F_full.segment<3>(e(0) * 4) * 0.5 + F_full.segment<3>(e(1) * 4) * 0.5;

    F.segment<3>(eidx * 3) += fu;
  });

  threadutils::for_each(0, num_faces, [&](int fidx) {
    const int gidx = fidx + num_edges;
    const auto& f = m_faces.row(fidx);
    const Vector3s& fw = m_face_weights[fidx];

    Vector3s fu = F_full.segment<3>(f[0] * 4) * fw[0] +
                  F_full.segment<3>(f[1] * 4) * fw[1] +
                  F_full.segment<3>(f[2] * 4) * fw[2];

    F.segment<3>(gidx * 3) += fu;
  });
}

const std::vector<Vector3s>& TwoDScene::getFaceWeights() const {
  return m_face_weights;
}

void TwoDScene::accumulateFluidNodeGradU(std::vector<VectorXs>& node_rhs_x,
                                         std::vector<VectorXs>& node_rhs_y,
                                         std::vector<VectorXs>& node_rhs_z,
                                         scalar coeff) {
  for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i) {
    if (m_forces[i]->flag() & 1)
      m_forces[i]->addLiquidGradEToNode(*this, node_rhs_x, node_rhs_y,
                                        node_rhs_z, coeff);
  }
}

void TwoDScene::accumulateGradU(VectorXs& F, const VectorXs& dx,
                                const VectorXs& dv) {
  assert(dx.size() == dv.size());

  if (F.size() == 0) return;

  VectorXs combined_mass = m_m + m_fluid_m;

  // Accumulate all energy gradients
  if (dx.size() == 0)
    for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i) {
      if (m_forces[i]->flag() & 1)
        m_forces[i]->addGradEToTotal(m_x, m_v, combined_mass, m_volume_fraction,
                                     m_sim_info.lambda, F);
    }
  else {
    VectorXs ddx = m_x + dx;
    VectorXs ddv = m_v + dv;

    for (std::vector<Force*>::size_type i = 0; i < m_forces.size(); ++i) {
      if (m_forces[i]->flag() & 1)
        m_forces[i]->addGradEToTotal(ddx, ddv, m_m, m_volume_fraction,
                                     m_sim_info.lambda, F);
    }
  }
}

bool TwoDScene::useAMGPCGSolid() const {
  return m_sim_info.use_amgpcg_solid;
}

scalar TwoDScene::totalFluidVolumeSoftElasto() const {
  const int num_elasto = getNumElastoParticles();
  const int num_fluids = getNumFluidParticles();
  return m_fluid_vol.segment(num_elasto, num_fluids).sum();
}

bool TwoDScene::isGaussFixed(int i) const {
  const int num_edges = getNumEdges();

  bool is_edge = i < num_edges;

  bool aisfixed = (is_edge && (isFixed(m_edges(i, 0)) & 1) &&
                   (isFixed(m_edges(i, 1)) & 1)) ||
                  (!is_edge && (isFixed(m_faces(i - num_edges, 0)) & 1) &&
                   (isFixed(m_faces(i - num_edges, 1)) & 1) &&
                   (isFixed(m_faces(i - num_edges, 2)) & 1));

  return aisfixed;
}

/*!
 * calculate velocity gradient on elements
 */
void TwoDScene::accumulateGaussGradU(MatrixXs& F, const VectorXs& dx,
                                     const VectorXs& dv) {
  //    MatrixXs dFe;
  //    dFe.resize(m_Fe_gauss.rows(), m_Fe_gauss.cols());

  //    std::cout<<dFe.norm()<<std::endl;
  assert(!std::isnan(m_dFe_gauss.sum()));

  //    std::cout<<m_vol_gauss<<std::endl;
  //    std::cout<<"dFe: "<<dFe<<std::endl;

  const int num_gauss = getNumGausses();
  const int num_edges = getNumEdges();
  threadutils::for_each(0, num_gauss, [&](int i) {
    bool is_edge = i < num_edges;
    bool is_cloth = !is_edge && i < getNumEdges() + getNumFaces();

    scalar psi_coeff = 1.0;
    if (is_edge || is_cloth) {
      psi_coeff = pow(m_volume_fraction_gauss[i], m_sim_info.lambda);
    }

    const Vector3s& dFe3 = m_dFe_gauss.block<3, 1>(i * 3, 2);
    const Vector3sT& d3t = m_d_gauss.block<3, 1>(i * 3, 2).transpose();
    F.block<3, 3>(3 * i, 0) += psi_coeff * m_rest_vol_gauss(i) * dFe3 *
                               d3t;  // * m_D_inv_gauss.block<3,3>(i*3,0);

    //        std::cout<<"vol: "<<m_vol_gauss(i)<<std::endl;
    //        std::cout<<"dFe3: "<< dFe3 <<std::endl;
    //        std::cout<<"d3t: "<<d3t<<std::endl;
    //
    if (is_edge) {
      const Vector3s& dFe2 = m_dFe_gauss.block<3, 1>(i * 3, 1);
      //            std::cout << "dFe2: " << dFe2 << std::endl;

      const Vector3sT& d2t = m_d_gauss.block<3, 1>(i * 3, 1).transpose();
      F.block<3, 3>(3 * i, 0) += psi_coeff * m_rest_vol_gauss(i) * dFe2 *
                                 d2t;  // * m_D_inv_gauss.block<3,3>(i*3,0);
    }

    //    std::cout << "F: " << i << "\n" << F.block<3,3>(3*i, 0) << std::endl;
  });
}

/*!
 * accumulate Hessian matrix
 */
void TwoDScene::accumulateddUdxdx(TripletXs& A, const scalar& dt, int base_idx,
                                  const VectorXs& dx, const VectorXs& dv) {
  assert(dx.size() == dv.size());

  int num_hess = base_idx;
  const int num_force = m_forces.size();

  std::vector<int> offsets(num_force);

  for (int i = 0; i < num_force; ++i) {
    offsets[i] = num_hess;
    num_hess += m_forces[i]->numHessX();
  }

  if ((int)A.size() != num_hess) A.resize(num_hess);

  if (dx.size() == 0) {
    threadutils::for_each(0, num_force, [&](int i) {
      if (!m_forces[i]->parallelized())
        m_forces[i]->addHessXToTotal(m_x, m_v, m_m, m_volume_fraction,
                                     m_sim_info.lambda, A, offsets[i], dt);
    });

    for (int i = 0; i < num_force; ++i) {
      if (m_forces[i]->parallelized())
        m_forces[i]->addHessXToTotal(m_x, m_v, m_m, m_volume_fraction,
                                     m_sim_info.lambda, A, offsets[i], dt);
    }
  } else {
    VectorXs idx = m_x + dx;
    VectorXs idv = m_v + dv;

    threadutils::for_each(0, num_force, [&](int i) {
      if (!m_forces[i]->parallelized())
        m_forces[i]->addHessXToTotal(idx, idv, m_m, m_volume_fraction,
                                     m_sim_info.lambda, A, offsets[i], dt);
    });

    for (int i = 0; i < num_force; ++i) {
      if (m_forces[i]->parallelized())
        m_forces[i]->addHessXToTotal(idx, idv, m_m, m_volume_fraction,
                                     m_sim_info.lambda, A, offsets[i], dt);
    }
  }
}

void TwoDScene::updateMultipliers(const scalar& dt) {
  const int num_force = m_forces.size();
  for (int i = 0; i < num_force; ++i) {
    m_forces[i]->updateMultipliers(m_x, m_v, m_m, m_volume_fraction,
                                   m_sim_info.lambda, dt);
  }
}

/*!
 * accumulate Hessian for the twisting DOF
 */
void TwoDScene::accumulateAngularddUdxdx(TripletXs& A, const scalar& dt,
                                         int base_idx, const VectorXs& dx,
                                         const VectorXs& dv) {
  assert(dx.size() == dv.size());

  int num_hess = base_idx;
  const int num_force = m_forces.size();

  std::vector<int> offsets(num_force);

  for (int i = 0; i < num_force; ++i) {
    offsets[i] = num_hess;
    num_hess += m_forces[i]->numAngularHessX();
  }

  if ((int)A.size() != num_hess) A.resize(num_hess);

  if (dx.size() == 0)
    for (int i = 0; i < num_force; ++i) {
      m_forces[i]->addAngularHessXToTotal(m_x, m_v, m_m, m_volume_fraction,
                                          m_sim_info.lambda, A, offsets[i],
                                          dt);
    }
  else {
    VectorXs idx = m_x + dx;
    VectorXs idv = m_v + dv;

    for (int i = 0; i < num_force; ++i) {
      m_forces[i]->addAngularHessXToTotal(idx, idv, m_m, m_volume_fraction,
                                          m_sim_info.lambda, A, offsets[i],
                                          dt);
    }
  }
}

void TwoDScene::dump_geometry(std::string filename) {
  int s = getNumParticles();
  std::ofstream myfile;
  myfile.open(filename);
  myfile << s << std::endl;
  for (int i = 0; i < s; i++) {
    myfile << m_x[4 * i] << " " << m_x[4 * i + 1] << " " << m_x[4 * i + 2]
           << std::endl;
  }
  myfile.close();
}

void TwoDScene::stepScript(const scalar& dt, const scalar& current_time) {
  threadutils::for_each(0, (int)m_scripts.size(), [&](int i) {
    m_scripts[i]->stepScript(dt, current_time);
  });
}

/*!
 * update rigid body level set
 */
void TwoDScene::updateSolidPhi() {

  m_particle_buckets.for_each_bucket([&](int bucket_idx) {
    if (!m_bucket_activated[bucket_idx]) return;

    const int num_nodes = getNumNodes(bucket_idx);

    VectorXs& node_phi = m_node_solid_phi[bucket_idx];

    VectorXs& node_solid_vel_x = m_node_solid_vel_x[bucket_idx];
    VectorXs& node_solid_vel_y = m_node_solid_vel_y[bucket_idx];
    VectorXs& node_solid_vel_z = m_node_solid_vel_z[bucket_idx];

    for (int i = 0; i < num_nodes; ++i) {
      Vector3s vel;
      node_phi(i) =
          computePhiVel(getNodePosSolidPhi(bucket_idx, i), vel);

      computePhiVel(getNodePosX(bucket_idx, i), vel);
      node_solid_vel_x(i) = vel(0);

      computePhiVel(getNodePosY(bucket_idx, i), vel);
      node_solid_vel_y(i) = vel(1);

      computePhiVel(getNodePosZ(bucket_idx, i), vel);
      node_solid_vel_z(i) = vel(2);
    }
  });
}

bool TwoDScene::isBucketActivated(int bucket_index) const {
  return m_bucket_activated[bucket_index];
}

/*!
 * apply scripts to transform objects
 */
void TwoDScene::applyScript(const scalar& dt) {
  const int np = getNumParticles();
  threadutils::for_each(0, np, [&](int i) {
    if (!(isFixed(i))) return;

    const int sg_idx = m_particle_group[i];
    if (sg_idx < 0 || sg_idx >= (int)m_group_pos.size()) return;

    const Eigen::Quaternion<scalar>& q = m_group_rot[sg_idx];
    const Vector3s& t = m_group_pos[sg_idx];

    const Eigen::Quaternion<scalar>& q_prev = m_group_prev_rot[sg_idx];
    const Vector3s& t_prev = m_group_prev_pos[sg_idx];

    const Eigen::Quaternion<scalar> q_diff = q * q_prev.inverse();

    if (isFixed(i) & 1) {
      Eigen::Quaternion<scalar> p;
      const Vector3s x0 = m_rest_x.segment<3>(i * 4) - t_prev;

      Vector3s trans_x0 = q_diff * x0 + t;

      m_rest_x.segment<3>(i * 4) = trans_x0;

      if (m_particle_to_surfel[i] >= 0) {
        m_v.segment<3>(i * 4) = (trans_x0 - m_x.segment<3>(i * 4)) / dt;
      }
    }

    if (m_twist[i] && isFixed(i) & 2) {
      Vector3s dir = getRestTwistDir(i);
      m_rest_x(i * 4 + 3) += mathutils::twist_component(q_diff, dir);
    }
  });

  const int num_surfels = getNumSurfels();

  threadutils::for_each(0, num_surfels, [&](int i) {
    const int pidx = m_surfels[i];
    const int sg_idx = m_particle_group[pidx];

    const Eigen::Quaternion<scalar>& q = m_group_rot[sg_idx];

    const Eigen::Quaternion<scalar>& q_prev = m_group_prev_rot[sg_idx];

    const Eigen::Quaternion<scalar> q_diff = q * q_prev.inverse();

    m_surfel_norms[i] = q_diff * m_surfel_norms[i];
  });

}

void TwoDScene::checkConsistency() {
  // TODO: Add more checks
}
