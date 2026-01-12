/**
 * @file test_para.cpp
 * @brief Test the parameter parser.
 */

#define DEBUG 1
#include "../include/myprompt.hpp"
#include "../include/para.hpp"
using namespace std;
using namespace otf;

int main()
{
    runtime_para para("../examples/galotfa.toml");

    INFO("GLOBAL PARAMETERS:");
    if (para.enableOtf)
    {
        INFO("On the fly analysis is enabled.");
    }
    else
    {
        INFO("On the fly analysis is forbidden.");
    }

    INFO("Output to  [%s]/[%s]", para.outputDir.c_str(), para.fileName.c_str());
    INFO("Max iteration [%u], epsilon [%g]", para.maxIter, para.epsilon);

    INFO("ORBITAL LOG PARAMETERS:");
    if (para.orbit->enable())
    {
        INFO("Orbital log is enabled.");
    }
    INFO("Log period: %d", para.orbit->period());
    INFO("Particle types to be logged:");
    for (const auto& id : para.orbit->sampleTypes())
    {
        INFO("%d ", id);
    }

    INFO("Method of id determination: ");
    switch (para.orbit->method())
    {
    case orbit::id_selection_method::RANDOM:
        INFO("Random selection.");
        break;
    case orbit::id_selection_method::TXTFILE:
        INFO("By a text file of id list.");
        break;
    default:
        ERROR("Get into an unexpected branch!");
    }
    INFO("Random selection fraction: %g.", para.orbit->fraction());
    INFO("ID list filename : %s", para.orbit->idfile().c_str());

    if (para.orbit->recenter.enable)
    {
        INFO("Enable recenter before log of orbits.");
    }
    switch (para.orbit->recenter.method)
    {
    case otf::recenter_method::COM:
        INFO("By conter of mass.");
        break;
    case otf::recenter_method::MBP:
        INFO("By most bound particle.");
        break;
    default:
        ERROR("Get into an unexpected branch!");
    }
    INFO("Recenter radius: %g", para.orbit->recenter.radius);
    INFO("Initial guess of the recenter:");
    for (auto& coord : para.orbit->recenter.initialGuess)
    {
        INFO("%g ", coord);
    }
    INFO("Particle types to be used as recenter anchors:");
    for (auto& id : para.orbit->recenter.anchorIds)
    {
        INFO("%d ", id);
    }

    INFO("COMPONENT PARAMETERS");
    for (auto& comp : para.comps)
    {
        INFO("Get an component: %s", comp.first.c_str());
        INFO("Analysis period of this component: %d", comp.second->period);
        INFO("Particle types in this component:");
        for (auto& id : comp.second->types)
        {
            INFO("%d ", id);
        }
        if (comp.second->recenter.enable)
        {
            INFO("Recenter of this component is enabled.");
        }
        string method = "";
        switch (comp.second->recenter.method)
        {
        case otf::recenter_method::COM:
            method = "center of mass";
            break;
        case otf::recenter_method::MBP:
            method = "most bound particle";
            break;
        default:
            ERROR("Get into an unexpected branch!");
        }
        INFO("Potential method for recenter: %s.", method.c_str());
        INFO("Potential enclosed radius for recenter: %g.",
             comp.second->recenter.radius);
        INFO("Initial guess of the recenter:");
        for (auto& coord : para.orbit->recenter.initialGuess)
        {
            INFO("%g ", coord);
        }
        if (comp.second->align.enable)
        {
            INFO("Alignment of this component is enabled.");
            INFO("Potential enclosed radius for align: %g.",
                 comp.second->align.radius);
        }
        if (comp.second->image.enable)
        {
            INFO("Image of this component is enabled.");
            INFO("Image half length: %g.", comp.second->image.halfLength);
            INFO("Image bin number: %d.", comp.second->image.binNum);
        }
        if (comp.second->sBar.enable)
        {
            INFO("A2 of this component is enabled.");
            INFO("A2 rmin : %g.", comp.second->sBar.rmin);
            INFO("A2 rmax : %g.", comp.second->sBar.rmax);
        }

        if (comp.second->barAngle.enable)
        {
            INFO("barAngle of this component is enabled.");
            INFO("barAngle rmin : %g.", comp.second->barAngle.rmin);
            INFO("barAngle rmax : %g.", comp.second->barAngle.rmax);
        }

        if (comp.second->sBuckle.enable)
        {
            INFO("buckle of this component is enabled.");
            INFO("buckle rmin : %g.", comp.second->sBuckle.rmin);
            INFO("buckle rmax : %g.", comp.second->sBuckle.rmax);
        }

        if (comp.second->A2profile.enable)
        {
            INFO("Radial A2 profile of this component is enabled.");
            INFO("Radial A2 profile rmin : %g.", comp.second->A2profile.rmin);
            INFO("Radial A2 profile rmax : %g.", comp.second->A2profile.rmax);
            INFO("Radial A2 profile binnum : %u.",
                 comp.second->A2profile.binNum);
        }
    }

    return 0;
}
