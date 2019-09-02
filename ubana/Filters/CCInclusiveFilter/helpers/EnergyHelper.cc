#include "EnergyHelper.h"

EnergyHelper::EnergyHelper()
{
}

void EnergyHelper::reconfigure(fhicl::ParameterSet const &pset)
{
    _recombination_factor = pset.get<float>("RecombinationFactor", 0.62);
    _dQdx_rectangle_length = pset.get<float>("dQdxRectangleLength", 4);
    _dQdx_rectangle_width = pset.get<float>("dQdxRectangleWidth", 1);
    _betap = pset.get<float>("RecombinationBeta", 0.212);
    _alpha = pset.get<float>("RecombinationAlpha", 0.93);
    m_isData = pset.get<bool>("is_data", false);

    if(m_isData){
        _gain = _data_gain;
    }
    else
    {
        _gain=_mc_gain;
    }
}

void EnergyHelper::energy_from_hits(const lar_pandora::ClusterVector &clusters,
                                    std::vector<uint> &nHits,
                                    std::vector<float> &pfenergy)
{

    nHits.resize(3);
    pfenergy.resize(3);
    for (auto cluster : clusters)
    {
        if (cluster->isValid())
        {
            // https://arxiv.org/pdf/1704.02927.pdf
            nHits[cluster->View()]= cluster->NHits();
            pfenergy[cluster->View()] = cluster->Integral() * _gain[cluster->View()] * _work_function / _recombination_factor / 1000; // convert MeV to GeV
        }
    }
}

void EnergyHelper::dQdx(const TVector3 &pfp_dir,
                        const lar_pandora::ClusterVector &clusters,
                        const lar_pandora::ClustersToHits &hits_per_cluster,
                        std::vector<float> &dqdx,
                        std::vector<std::vector<float>> &dqdx_hits,
                        std::vector<float> &pitches)
{

    float tolerance = 0.001;

    for (auto _cl : clusters)
    {
        std::vector<art::Ptr<recob::Hit>> hits = hits_per_cluster.at(_cl);

        std::vector<float> cluster_axis;
        std::vector<float> cluster_start;
        std::vector<float> cluster_end;
        float start_x = _detprop->ConvertTicksToX(_cl->StartTick(), _cl->Plane());
        float end_x = _detprop->ConvertTicksToX(_cl->EndTick(), _cl->Plane());
        float pitch = getPitch(pfp_dir, _cl->Plane().Plane);
        if (pitch >= 0)
        {
            std::reverse(hits.begin(), hits.end());
            cluster_axis = {cos(_cl->StartAngle()),
                            sin(_cl->StartAngle())};

            cluster_start = {_cl->StartWire() * _wire_spacing - tolerance * cos(_cl->StartAngle()),
                             start_x - tolerance * sin(_cl->StartAngle())};
            cluster_end = {_cl->EndWire() * _wire_spacing, end_x};
        }
        else
        {
            cluster_axis = {-1 * cos(_cl->StartAngle()),
                            -1 * sin(_cl->StartAngle())};
            cluster_start = {_cl->EndWire() * _wire_spacing + tolerance * cos(_cl->StartAngle()),
                             end_x + tolerance * sin(_cl->StartAngle())};
            cluster_end = {_cl->StartWire() * _wire_spacing, start_x};
        }

        float cluster_length = sqrt(pow(cluster_end[0] - cluster_start[0], 2) +
                                    pow(cluster_end[1] - cluster_start[1], 2));
        if (cluster_length <= 0)
            continue;

        // Build rectangle 4 x 1 cm around the cluster axis
        std::vector<std::vector<float>> points = buildRectangle(_dQdx_rectangle_length, _dQdx_rectangle_width, cluster_start, cluster_axis);

        std::vector<float> dqdxs;
        for (auto &hit : hits)
        {
            float hit_pos_x = _detprop->ConvertTicksToX(hit->PeakTime(), _cl->Plane());
            std::vector<float> hit_pos = {hit->WireID().Wire * _wire_spacing, hit_pos_x};

            bool is_within = isInside(hit_pos, points);
            if (is_within)
            {
                float q = hit->Integral() * _gain[_cl->Plane().Plane];
                dqdxs.push_back(q / fabs(pitch));
                dqdx_hits[_cl->Plane().Plane].push_back(q / fabs(pitch));
            }
        }

        // Get the median
        if (dqdxs.size() > 0)
        {
            std::nth_element(dqdxs.begin(), dqdxs.begin() + dqdxs.size() / 2, dqdxs.end());
            dqdx[_cl->Plane().Plane] = dqdxs[dqdxs.size() / 2];
            pitches[_cl->Plane().Plane] = pitch;
        }
    }
}

float EnergyHelper::getPitch(const TVector3 &direction,
                             const int &pl) const
{
    // prepare a direction vector for the plane
    TVector3 wireDir = {0., 0., 0.};
    // the direction of the plane is the vector uniting two consecutive wires
    // such that this vector is perpendicular to both wires
    // basically this is the vector perpendicular to the wire length direction,
    // and still in the wire-plane direction
    if (pl == 0)
        wireDir = {0., -sqrt(3) / 2., 1 / 2.};
    else if (pl == 1)
        wireDir = {0., sqrt(3) / 2., 1 / 2.};
    else if (pl == 2)
        wireDir = {0., 0., 1.};

    // cosine between shower direction and plane direction gives the factor
    // by which to divide 0.3, the minimum wire-spacing
    float minWireSpacing = 0.3;
    float cos = wireDir.Dot(direction);

    cos /= (wireDir.Mag() * direction.Mag());
    // if cosine is 0 the direction is perpendicular and the wire-spacing is
    // infinite
    if (cos == 0)
        return std::numeric_limits<float>::max();
    float pitch = minWireSpacing / cos;
    return pitch;
}

std::vector<std::vector<float>> EnergyHelper::buildRectangle(float length, float width,
                                                             std::vector<float> &start,
                                                             std::vector<float> &axis)
{
    float perp_axis[2] = {-axis[1], axis[0]};

    std::vector<float> p1 = {start[0] + perp_axis[0] * width / 2,
                             start[1] + perp_axis[1] * width / 2};

    std::vector<float> p2 = {p1[0] + axis[0] * length,
                             p1[1] + axis[1] * length};

    std::vector<float> p3 = {start[0] - perp_axis[0] * width / 2,
                             start[1] - perp_axis[1] * width / 2};

    std::vector<float> p4 = {p3[0] + axis[0] * length,
                             p3[1] + axis[1] * length};
    std::vector<std::vector<float>> points = {p1, p2, p4, p3};
    return points;
}

bool EnergyHelper::isInside(std::vector<float> P,
                            std::vector<std::vector<float>> V)
{

    int nvert = (int)V.size();

    int i, j, c = 0;
    for (i = 0, j = nvert - 1; i < nvert; j = i++)
    {
        if (((V[i][1] > P[1]) != (V[j][1] > P[1])) &&
            (P[0] < (V[j][0] - V[i][0]) * (P[1] - V[i][1]) / (V[j][1] - V[i][1]) + V[i][0]))
            c = !c;
    }
    return c;
}

std::vector<float> EnergyHelper::dEdx_from_dQdx(std::vector<float> dqdx)
{
    std::vector<float> dedx;

    for (size_t i = 0; i < dqdx.size(); i++)
    {
        // More advance model
        // float Rho = 1.4;
        // float Efield = 0.273;
        // dedx.push_back((exp(dqdx[i] * (_betap / (Rho * Efield)) * _work_function) - _alpha) / (_betap / (Rho * Efield)));
        
        // Fake linear model: https://arxiv.org/pdf/1704.02927.pdf pag. 16
        dedx.push_back( dqdx[i] * _work_function / _recombination_factor);
    }
    return dedx;
}