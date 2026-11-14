#ifndef XZZ_RJR_RECONSTRUCTION_H
#define XZZ_RJR_RECONSTRUCTION_H

#include "RestFrames/RestFrames.hh"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <cmath>
#include <stdexcept>

namespace XZZRJR {

struct Result {
    double mass = -999.;
    double visibleMass = -999.;
    double invisibleMass = -999.;
    // 0: success; 1: invalid input; 2: reconstruction failed; 3: invalid output.
    int status = 1;
};

// One instance per worker/thread; never share an instance concurrently.
class Reconstruction {
public:
    Reconstruction()
        : lab_("RJR_LAB", "LAB"), x_("RJR_X", "X"),
          za_("RJR_Za", "Za"), zb_("RJR_Zb", "Zb"),
          l1_("RJR_L1", "L1"), l2_("RJR_L2", "L2"),
          n1_("RJR_N1", "N1"), n2_("RJR_N2", "N2"),
          inv_("RJR_INV", "INV"), contra_("RJR_Contra", "Contra") {
        lab_.SetChildFrame(x_);
        x_.AddChildFrame(za_);
        x_.AddChildFrame(zb_);
        za_.AddChildFrame(l1_);
        za_.AddChildFrame(l2_);
        zb_.AddChildFrame(n1_);
        zb_.AddChildFrame(n2_);
        if (!lab_.InitializeTree())
            throw std::runtime_error("RJR: InitializeTree failed");
        inv_.AddFrame(n1_);
        inv_.AddFrame(n2_);
        // Match the contra-boost-only reference, with no rapidity or fixed-mass rule.
        inv_.AddJigsaw(contra_);
        contra_.AddVisibleFrame(l1_, 0);
        contra_.AddVisibleFrame(l2_, 1);
        contra_.AddInvisibleFrame(n1_, 0);
        contra_.AddInvisibleFrame(n2_, 1);
        if (!lab_.InitializeAnalysis())
            throw std::runtime_error("RJR: InitializeAnalysis failed");
    }

    Reconstruction(const Reconstruction&) = delete;
    Reconstruction& operator=(const Reconstruction&) = delete;
    Reconstruction(Reconstruction&&) = delete;
    Reconstruction& operator=(Reconstruction&&) = delete;

    Result reconstruct(const TLorentzVector& a, const TLorentzVector& b,
                       double metPt, double metPhi) {
        Result result;
        lab_.ClearEvent();
        if (!valid(a) || !valid(b) || !std::isfinite(metPt) ||
            !std::isfinite(metPhi) || metPt < 0.)
            return result;
        l1_.SetLabFrameFourVector(a);
        l2_.SetLabFrameFourVector(b);
        inv_.SetLabFrameThreeVector(
            TVector3(metPt * std::cos(metPhi), metPt * std::sin(metPhi), 0.));
        if (!lab_.AnalyzeEvent()) {
            result.status = 2;
            return result;
        }
        const double mass = x_.GetMass();
        const double visible = za_.GetMass();
        const double invisible = zb_.GetMass();
        if (!std::isfinite(mass) || !std::isfinite(visible) ||
            !std::isfinite(invisible) || mass <= 0. ||
            visible < 0. || invisible < 0.) {
            result.status = 3;
            return result;
        }
        result.mass = mass;
        result.visibleMass = visible;
        result.invisibleMass = invisible;
        result.status = 0;
        return result;
    }

private:
    static bool valid(const TLorentzVector& p) {
        return std::isfinite(p.Px()) && std::isfinite(p.Py()) &&
               std::isfinite(p.Pz()) && std::isfinite(p.E()) &&
               std::isfinite(p.M2()) && p.E() > 0. && p.M2() >= 0.;
    }

    RestFrames::LabRecoFrame lab_;
    RestFrames::DecayRecoFrame x_, za_, zb_;
    RestFrames::VisibleRecoFrame l1_, l2_;
    RestFrames::InvisibleRecoFrame n1_, n2_;
    RestFrames::InvisibleGroup inv_;
    RestFrames::ContraBoostInvJigsaw contra_;
};

} // namespace XZZRJR
#endif

