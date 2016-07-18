// -*- C++ -*-
// KimParameterProvider.cpp:  Get parameters from a string (for OpenKIM).
//
// Copyright (C) 2012-2013 Jakob Schiotz and the Department of Physics,
// Technical University of Denmark.  Email: schiotz@fysik.dtu.dk
//
// This file is part of Asap version 3.
// Asap is released under the GNU Lesser Public License (LGPL) version 3.
// However, the parts of Asap distributed within the OpenKIM project
// (including this file) are also released under the Common Development
// and Distribution License (CDDL) version 1.0.
//
// This program is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// version 3 as published by the Free Software Foundation.  Permission
// to use other versions of the GNU Lesser General Public License may
// granted by Jakob Schiotz or the head of department of the
// Department of Physics, Technical University of Denmark, as
// described in section 14 of the GNU General Public License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// and the GNU Lesser Public License along with this program.  If not,
// see <http://www.gnu.org/licenses/>.


#include "KimParameterProvider.h"
#include "Asap.h"
#include "Exception.h"
// #define ASAPDEBUG
#include "Debug.h"
#include "KIM_API_C.h"
#include "KIM_API_status.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <set>
#include <stdlib.h>

KimParameterProvider::KimParameterProvider(const char *parameter_filename,
                                           intptr_t *pkim)
{
  DEBUGPRINT;
  // Get a list of elements in the simulation.
  int nelements;
  int strlen;
  int ier;
  ier = KIM_API_get_num_sim_species(pkim, &nelements, &strlen);
  if (ier != KIM_STATUS_OK || nelements < 1)
    {
      std::cerr << "KIM_API_get_num_sim_species failed. " << ier << " " << nelements << std::endl;
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_num_sim_species failed.", ier);
      throw AsapError("Failed to get element list from OpenKIM.");
    }
  std::set<string> elements;
  for (int i = 0; i < nelements; i++)
    {
      const char* element;
      ier = KIM_API_get_sim_species(pkim, i, &element);
      string el = string(element);
      elements.insert(el);
      // std::cerr << "Simulation needs element *" << el << "*" << std::endl;
    }

  // Prepare to convert all parameters to other units.
  double energy_conv =  KIM_API_convert_to_act_unit(pkim, "A", "eV", "e",
      "K", "ps", 0.0, 1.0, 0.0, 0.0, 0.0, &ier);
  assert(ier == KIM_STATUS_OK);  // Should not be able to fail.
  double len_conv = KIM_API_convert_to_act_unit(pkim, "A", "eV", "e",
      "K", "ps", 1.0, 0.0, 0.0, 0.0, 0.0, &ier);
  assert(ier == KIM_STATUS_OK);  // Should not be able to fail.
  double inv_len_conv = KIM_API_convert_to_act_unit(pkim, "A", "eV", "e",
      "K", "ps", -1.0, 0.0, 0.0, 0.0, 0.0, &ier);
  assert(ier == KIM_STATUS_OK);  // Should not be able to fail.
  double inv_vol_conv = KIM_API_convert_to_act_unit(pkim, "A", "eV", "e",
      "K", "ps", -3.0, 0.0, 0.0, 0.0, 0.0, &ier);
  assert(ier == KIM_STATUS_OK);  // Should not be able to fail.

  // Read the parameter file.
  std::fstream pstr(parameter_filename, std::fstream::in);
  int numread = -1;
  emt_parameters *p = NULL;
  while(true)
    {
      string word;
      double value;
      pstr >> word >> value;
      //std::cerr << "read(" << word << ")(" << value << ")" << std::endl;
      numread++;
      if (word == "NEWELEMENT")
        {
          string name;
          pstr >> name;
          numread = 0;
          p = new emt_parameters;
          p->Z = (int) value;
          p->gamma1 = 0.0;        // These are calculated later!
          p->gamma2 = 0.0;
          p->name = name;
        }
      else if (word == "V0")
        p->V0 = value * energy_conv;
      else if (word == "kappa")
        p->kappa = value * inv_len_conv;
      else if (word == "eta2")
        p->eta2 = value * inv_len_conv;
      else if (word == "n0")
        p->neq = value * inv_vol_conv;
      else if (word == "S0")
        p->seq = value * len_conv;
      else if (word == "E0")
        p->e0 = value * energy_conv;
      else if (word == "lambda")
        p->lambda = value * inv_len_conv;
      else if (word == "mass")
        {
          p->mass = value;  // Never really used.
          p->invmass = 1.0/value;
        }
      else if (word == "ENDELEMENT")
        {
          if (numread != 9)
            {
              std::cerr << "Wrong number of parameters for element " << p->name << std::endl;
              throw AsapError("Error in parameterfile");
            }
          if (elements.count(p->name))
            {
              // We need this element
              p->index = params.size();  // Count the elements, starting at zero.
              params.push_back(p);
            }
          else
            delete p;  // We don't need it.
          p = NULL;
        }
      else if (word == "STOP")
        break;
      else
        {
          std::cerr << "Unknown keyword in parameter file: " << word << std::endl;
          throw AsapError("Error in parameterfile");
        }
    }
}

emt_parameters *KimParameterProvider::GetNewParameters(int element)
{
  throw AsapError("Element missing from parameter file or unannounced by simulation: number = ") << element;
}

