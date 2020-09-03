//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "SpringForce.h"

SpringForce::SpringForce(const Vector2iT& endpoints, const scalar& k,
                         const scalar& l0, const scalar& b)
    : Force(),
      m_endpoints(endpoints),
      m_k(k),
      m_l0(l0),
      m_b(b),
      m_zero_length(l0 == 0.0) {
  assert(m_k >= 0.0);
  assert(m_l0 >= 0.0);
  assert(m_b >= 0.0);
}

SpringForce::~SpringForce() {}

void SpringForce::addEnergyToTotal(const VectorXs& x, const VectorXs& v,
                                   const VectorXs& m, const VectorXs& psi,
                                   const scalar& lambda, scalar& E) {
  assert(x.size() == v.size());
  assert(x.size() % 4 == 0);

  scalar l =
      (x.segment<3>(4 * m_endpoints(1)) - x.segment<3>(4 * m_endpoints(0)))
          .norm();
  scalar psi_coeff =
      pow((psi(m_endpoints(1)) + psi(m_endpoints(0))) * 0.5, lambda);
  E += 0.5 * m_k * (l - m_l0) * (l - m_l0) * psi_coeff;
}

void SpringForce::addGradEToTotal(const VectorXs& x, const VectorXs& v,
                                  const VectorXs& m, const VectorXs& psi,
                                  const scalar& lambda, VectorXs& gradE) {
  assert(x.size() == v.size());
  assert(x.size() == gradE.size());
  assert(x.size() % 4 == 0);
  scalar psi_coeff =
      pow((psi(m_endpoints(1)) + psi(m_endpoints(0))) * 0.5, lambda);

  if (m_zero_length) {
    Vector3s nhat =
        x.segment<3>(4 * m_endpoints(1)) - x.segment<3>(4 * m_endpoints(0));
    nhat *= m_k * psi_coeff;
    gradE.segment<3>(4 * m_endpoints(0)) -= nhat;
    gradE.segment<3>(4 * m_endpoints(1)) += nhat;
  } else {
    // Compute the elastic component
    Vector3s nhat =
        x.segment<3>(4 * m_endpoints(1)) - x.segment<3>(4 * m_endpoints(0));
    scalar l = nhat.norm();
    assert(l != 0.0);
    nhat /= l;
    nhat *= m_k * (l - m_l0) * psi_coeff;
    gradE.segment<3>(4 * m_endpoints(0)) -= nhat;
    gradE.segment<3>(4 * m_endpoints(1)) += nhat;
  }
}

void SpringForce::addHessXToTotal(const VectorXs& x, const VectorXs& v,
                                  const VectorXs& m, const VectorXs& psi,
                                  const scalar& lambda, TripletXs& hessE,
                                  int hessE_index, const scalar& dt) {
  assert(x.size() == v.size());
  assert(x.size() == m.size());
  assert(x.size() % 4 == 0);
  scalar psi_coeff =
      pow((psi(m_endpoints(1)) + psi(m_endpoints(0))) * 0.5, lambda);

  // Contribution from elastic component
  if (m_zero_length) {
    const scalar coeff = m_k * psi_coeff;
    for (int r = 0; r < 3; ++r) {
      hessE[hessE_index + r * 4 + 0] =
          Triplets(4 * m_endpoints(0) + r, 4 * m_endpoints(0) + r, coeff);
      hessE[hessE_index + r * 4 + 1] =
          Triplets(4 * m_endpoints(1) + r, 4 * m_endpoints(1) + r, coeff);
      hessE[hessE_index + r * 4 + 2] =
          Triplets(4 * m_endpoints(0) + r, 4 * m_endpoints(1) + r, -coeff);
      hessE[hessE_index + r * 4 + 3] =
          Triplets(4 * m_endpoints(1) + r, 4 * m_endpoints(0) + r, -coeff);
    }
  } else {
    Vector3s nhat =
        x.segment<3>(4 * m_endpoints(1)) - x.segment<3>(4 * m_endpoints(0));
    scalar l = nhat.norm();
    assert(l != 0);
    nhat /= l;

    Matrix3s hess;
    hess = nhat * nhat.transpose();
    hess += (l - m_l0) * (Matrix3s::Identity() - hess) / l;
    hess *= m_k * psi_coeff;

    for (int r = 0; r < 3; ++r)
      for (int s = 0; s < 3; ++s) {
        int element_idx = r * 3 + s;
        hessE[hessE_index + element_idx * 4 + 0] = Triplets(
            4 * m_endpoints(0) + r, 4 * m_endpoints(0) + s, hess(r, s));
        hessE[hessE_index + element_idx * 4 + 1] = Triplets(
            4 * m_endpoints(1) + r, 4 * m_endpoints(1) + s, hess(r, s));
        hessE[hessE_index + element_idx * 4 + 2] = Triplets(
            4 * m_endpoints(0) + r, 4 * m_endpoints(1) + s, -hess(r, s));
        hessE[hessE_index + element_idx * 4 + 3] = Triplets(
            4 * m_endpoints(1) + r, 4 * m_endpoints(0) + s, -hess(r, s));
      }
  }
}

void SpringForce::updateMultipliers(const VectorXs& x, const VectorXs& vplus,
                                    const VectorXs& m, const VectorXs& psi,
                                    const scalar& lambda, const scalar& dt) {}

void SpringForce::preCompute() {}

void SpringForce::updateStartState() {}

int SpringForce::numHessX() {
  if (m_zero_length) {
    return 3 * 4;
  } else {
    return 9 * 4;
  }
}

Force* SpringForce::createNewCopy() { return new SpringForce(*this); }

int SpringForce::flag() const { return 1; }
