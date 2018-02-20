/*
 *  Player - One Hell of a Robot Server
 *  Copyright (C) 2000  Brian Gerkey   &  Kasper Stoy
 *                      gerkey@usc.edu    kaspers@robotics.usc.edu
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
///////////////////////////////////////////////////////////////////////////
//
// Desc: AMCL odometry routines
// Author: Andrew Howard
// Date: 6 Feb 2003
// CVS: $Id: amcl_odom.cc 7057 2008-10-02 00:44:06Z gbiggs $
//
///////////////////////////////////////////////////////////////////////////

#include <algorithm>

#include <sys/types.h> // required by Darwin
#include <math.h>
#include <assert.h>
#include "amcl/sensors/amcl_odom.h"

using namespace amcl;

static double
normalize(double z)
{
  return atan2(sin(z),cos(z));
}
static double
angle_diff(double a, double b)
{
  double d1, d2;
  a = normalize(a);
  b = normalize(b);
  d1 = a-b;
  d2 = 2*M_PI - fabs(d1);
  if(d1 > 0)
    d2 *= -1.0;
  if(fabs(d1) < fabs(d2))
    return(d1);
  else
    return(d2);
}

////////////////////////////////////////////////////////////////////////////////
// Default constructor
AMCLOdom::AMCLOdom() : AMCLSensor()
{
  this->time = 0.0;
}

void
AMCLOdom::SetModelDiff(double alpha1, 
                       double alpha2, 
                       double alpha3, 
                       double alpha4)
{
  this->model_type = ODOM_MODEL_DIFF;
  this->alpha1 = alpha1;
  this->alpha2 = alpha2;
  this->alpha3 = alpha3;
  this->alpha4 = alpha4;
}

void
AMCLOdom::SetModelOmni(double alpha1, 
                       double alpha2, 
                       double alpha3, 
                       double alpha4,
                       double alpha5)
{
  this->model_type = ODOM_MODEL_OMNI;
  this->alpha1 = alpha1;
  this->alpha2 = alpha2;
  this->alpha3 = alpha3;
  this->alpha4 = alpha4;
  this->alpha5 = alpha5;
}

void
AMCLOdom::SetModel( odom_model_t type,
                    double alpha1,
                    double alpha2,
                    double alpha3,
                    double alpha4,
                    double alpha5 )
{
  this->model_type = type;
  this->alpha1 = alpha1;
  this->alpha2 = alpha2;
  this->alpha3 = alpha3;
  this->alpha4 = alpha4;
  this->alpha5 = alpha5;
}

////////////////////////////////////////////////////////////////////////////////
// Apply the odometry sensor model for histogram of Grid Localization
double AMCLOdom::UpdateSensor(pf_t *grid, AMCLSensorData *data)
{
  double total = 0.0;
  AMCLOdomData* ndata;
  ndata = (AMCLOdomData*) data;
  pf_sample_set_t* set_a;
  pf_sample_set_t* set_b;
  pf_sample_t* sample_a;
  pf_sample_t* sample_b;
  pf_vector_t pose, delta;
  pf_vector_t old_pose = pf_vector_sub(ndata->pose, ndata->delta);
  set_a = grid->sets + grid->current_set;
  set_b = grid->sets + (grid->current_set + 1)%2;
  //calculate delta of ndata
  double delta_rot1, delta_trans, delta_rot2;
  double delta_rot1_hat, delta_trans_hat, delta_rot2_hat;
  if(sqrt(ndata->delta.v[1]*ndata->delta.v[1] + 
          ndata->delta.v[0]*ndata->delta.v[0]) < 0.01)
    delta_rot1 = 0.0;
  else
    delta_rot1 = angle_diff(atan2(ndata->delta.v[1], ndata->delta.v[0]),
                            old_pose.v[2]);
  delta_trans = sqrt(ndata->delta.v[0]*ndata->delta.v[0] +
                     ndata->delta.v[1]*ndata->delta.v[1]);
  delta_rot2 = angle_diff(ndata->delta.v[2], delta_rot1);
  assert(isnan(delta_rot1)==false);
  assert(isnan(delta_trans)==false);
  assert(isnan(delta_rot2)==false);
  //TODO parallelization
  //TODO struct
  bool screening = false;
  double min_trans, max_trans;
  double min_rot1, max_rot1;
  double epson = 1.0/set_a->sample_count;
  for (int i = 0 ; i < set_a->sample_count; ++i)
  {
    sample_a = set_a->samples + i;
    //TODO slective update
    if ( set_b->samples[i].weight < epson )//if bel(xt-1)
    {
      sample_a->weight = epson;//then bel(xt)~=epson
      total +=sample_a->weight;
      continue;
    }
    double probi = 0.0;
    double pmotion = 0.0;
    double a1,a2,a3;
    double b1,b2,b3;
    for (int j = 0 ; j < set_b->sample_count; ++j)
    {
      sample_b = set_b->samples + j;
      if ( sample_b->weight < epson )//if bel(xt-1)
      {
        continue;
      }
      delta.v[0] = sample_a->pose.v[0]-sample_b->pose.v[0];
      delta.v[1] = sample_a->pose.v[1]-sample_b->pose.v[1];
      //intput sample_a sample_b ndata, output weight
      //calculate delta of sample_a and sample_b
      delta_trans_hat = sqrt(delta.v[0]*delta.v[0] +
                         delta.v[1]*delta.v[1]);
      if(delta_trans_hat < 0.01)
        delta_rot1_hat = 0.0;
      else
        delta_rot1_hat = angle_diff(atan2(delta.v[1], delta.v[0]),
                                sample_b->pose.v[2]);

      a1 = delta_rot1 - delta_rot1_hat; 
      b1 = alpha1*delta_rot1_hat*delta_rot1_hat + 
          alpha2*delta_trans_hat*delta_trans_hat;
      //when a is larger than 4*sqrt(b), where sqrt(b) is standard deviation
      if(b1 != 0 && a1*a1>=16.0*b1)
        continue;

      delta.v[2] = angle_diff(sample_a->pose.v[2],sample_b->pose.v[2]);
      delta_rot2_hat = angle_diff(delta.v[2], delta_rot1_hat);

      a2 = delta_trans - delta_trans_hat; 
      b2 = alpha3*delta_trans_hat*delta_trans_hat + 
          alpha4*delta_rot1_hat*delta_rot1_hat +
          alpha4*delta_rot2_hat*delta_rot2_hat;
      a3 = delta_rot2 - delta_rot2_hat; 
      b3 = alpha1*delta_rot2_hat*delta_rot2_hat + 
          alpha2*delta_trans_hat*delta_trans_hat;
      if(a2*a2>=16*b2 || a3*a3>=16*b3 )
        continue;
      pmotion = pf_normal_distribution(a1,b1)*
                pf_normal_distribution(a2,b2)*
                pf_normal_distribution(a3,b3);
      probi += (sample_b->weight * pmotion);
      
    }
    sample_a->weight = probi;
    total +=sample_a->weight;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Apply the action model
bool AMCLOdom::UpdateAction(pf_t *pf, AMCLSensorData *data)
{
  AMCLOdomData *ndata;
  ndata = (AMCLOdomData*) data;

  // Compute the new sample poses
  pf_sample_set_t *set;

  set = pf->sets + pf->current_set;
  pf_vector_t old_pose = pf_vector_sub(ndata->pose, ndata->delta);

  switch( this->model_type )
  {
  case ODOM_MODEL_OMNI:
  {
    double delta_trans, delta_rot, delta_bearing;
    double delta_trans_hat, delta_rot_hat, delta_strafe_hat;

    delta_trans = sqrt(ndata->delta.v[0]*ndata->delta.v[0] +
                       ndata->delta.v[1]*ndata->delta.v[1]);
    delta_rot = ndata->delta.v[2];

    // Precompute a couple of things
    double trans_hat_stddev = (alpha3 * (delta_trans*delta_trans) +
                               alpha1 * (delta_rot*delta_rot));
    double rot_hat_stddev = (alpha4 * (delta_rot*delta_rot) +
                             alpha2 * (delta_trans*delta_trans));
    double strafe_hat_stddev = (alpha1 * (delta_rot*delta_rot) +
                                alpha5 * (delta_trans*delta_trans));

    for (int i = 0; i < set->sample_count; i++)
    {
      pf_sample_t* sample = set->samples + i;

      delta_bearing = angle_diff(atan2(ndata->delta.v[1], ndata->delta.v[0]),
                                 old_pose.v[2]) + sample->pose.v[2];
      double cs_bearing = cos(delta_bearing);
      double sn_bearing = sin(delta_bearing);

      // Sample pose differences
      delta_trans_hat = delta_trans + pf_ran_gaussian(trans_hat_stddev);
      delta_rot_hat = delta_rot + pf_ran_gaussian(rot_hat_stddev);
      delta_strafe_hat = 0 + pf_ran_gaussian(strafe_hat_stddev);
      // Apply sampled update to particle pose
      sample->pose.v[0] += (delta_trans_hat * cs_bearing + 
                            delta_strafe_hat * sn_bearing);
      sample->pose.v[1] += (delta_trans_hat * sn_bearing - 
                            delta_strafe_hat * cs_bearing);
      sample->pose.v[2] += delta_rot_hat ;
    }
  }
  break;
  case ODOM_MODEL_DIFF:
  {
    // Implement sample_motion_odometry (Prob Rob p 136)
    double delta_rot1, delta_trans, delta_rot2;
    double delta_rot1_hat, delta_trans_hat, delta_rot2_hat;
    double delta_rot1_noise, delta_rot2_noise;

    // Avoid computing a bearing from two poses that are extremely near each
    // other (happens on in-place rotation).
    if(sqrt(ndata->delta.v[1]*ndata->delta.v[1] + 
            ndata->delta.v[0]*ndata->delta.v[0]) < 0.01)
      delta_rot1 = 0.0;
    else
      delta_rot1 = angle_diff(atan2(ndata->delta.v[1], ndata->delta.v[0]),
                              old_pose.v[2]);
    delta_trans = sqrt(ndata->delta.v[0]*ndata->delta.v[0] +
                       ndata->delta.v[1]*ndata->delta.v[1]);
    delta_rot2 = angle_diff(ndata->delta.v[2], delta_rot1);

    // We want to treat backward and forward motion symmetrically for the
    // noise model to be applied below.  The standard model seems to assume
    // forward motion.
    delta_rot1_noise = std::min(fabs(angle_diff(delta_rot1,0.0)),
                                fabs(angle_diff(delta_rot1,M_PI)));
    delta_rot2_noise = std::min(fabs(angle_diff(delta_rot2,0.0)),
                                fabs(angle_diff(delta_rot2,M_PI)));

    for (int i = 0; i < set->sample_count; i++)
    {
      pf_sample_t* sample = set->samples + i;

      // Sample pose differences
      delta_rot1_hat = angle_diff(delta_rot1,
                                  pf_ran_gaussian(this->alpha1*delta_rot1_noise*delta_rot1_noise +
                                                  this->alpha2*delta_trans*delta_trans));
      delta_trans_hat = delta_trans - 
              pf_ran_gaussian(this->alpha3*delta_trans*delta_trans +
                              this->alpha4*delta_rot1_noise*delta_rot1_noise +
                              this->alpha4*delta_rot2_noise*delta_rot2_noise);
      delta_rot2_hat = angle_diff(delta_rot2,
                                  pf_ran_gaussian(this->alpha1*delta_rot2_noise*delta_rot2_noise +
                                                  this->alpha2*delta_trans*delta_trans));

      // Apply sampled update to particle pose
      sample->pose.v[0] += delta_trans_hat * 
              cos(sample->pose.v[2] + delta_rot1_hat);
      sample->pose.v[1] += delta_trans_hat * 
              sin(sample->pose.v[2] + delta_rot1_hat);
      sample->pose.v[2] += delta_rot1_hat + delta_rot2_hat;
    }
  }
  break;
  case ODOM_MODEL_OMNI_CORRECTED:
  {
    double delta_trans, delta_rot, delta_bearing;
    double delta_trans_hat, delta_rot_hat, delta_strafe_hat;

    delta_trans = sqrt(ndata->delta.v[0]*ndata->delta.v[0] +
                       ndata->delta.v[1]*ndata->delta.v[1]);
    delta_rot = ndata->delta.v[2];

    // Precompute a couple of things
    double trans_hat_stddev = sqrt( alpha3 * (delta_trans*delta_trans) +
                                    alpha1 * (delta_rot*delta_rot) );
    double rot_hat_stddev = sqrt( alpha4 * (delta_rot*delta_rot) +
                                  alpha2 * (delta_trans*delta_trans) );
    double strafe_hat_stddev = sqrt( alpha1 * (delta_rot*delta_rot) +
                                     alpha5 * (delta_trans*delta_trans) );

    for (int i = 0; i < set->sample_count; i++)
    {
      pf_sample_t* sample = set->samples + i;

      delta_bearing = angle_diff(atan2(ndata->delta.v[1], ndata->delta.v[0]),
                                 old_pose.v[2]) + sample->pose.v[2];
      double cs_bearing = cos(delta_bearing);
      double sn_bearing = sin(delta_bearing);

      // Sample pose differences
      delta_trans_hat = delta_trans + pf_ran_gaussian(trans_hat_stddev);
      delta_rot_hat = delta_rot + pf_ran_gaussian(rot_hat_stddev);
      delta_strafe_hat = 0 + pf_ran_gaussian(strafe_hat_stddev);
      // Apply sampled update to particle pose
      sample->pose.v[0] += (delta_trans_hat * cs_bearing + 
                            delta_strafe_hat * sn_bearing);
      sample->pose.v[1] += (delta_trans_hat * sn_bearing - 
                            delta_strafe_hat * cs_bearing);
      sample->pose.v[2] += delta_rot_hat ;
    }
  }
  break;
  case ODOM_MODEL_DIFF_CORRECTED:
  {
    // Implement sample_motion_odometry (Prob Rob p 136)
    double delta_rot1, delta_trans, delta_rot2;
    double delta_rot1_hat, delta_trans_hat, delta_rot2_hat;
    double delta_rot1_noise, delta_rot2_noise;

    // Avoid computing a bearing from two poses that are extremely near each
    // other (happens on in-place rotation).
    if(sqrt(ndata->delta.v[1]*ndata->delta.v[1] + 
            ndata->delta.v[0]*ndata->delta.v[0]) < 0.01)
      delta_rot1 = 0.0;
    else
      delta_rot1 = angle_diff(atan2(ndata->delta.v[1], ndata->delta.v[0]),
                              old_pose.v[2]);
    delta_trans = sqrt(ndata->delta.v[0]*ndata->delta.v[0] +
                       ndata->delta.v[1]*ndata->delta.v[1]);
    delta_rot2 = angle_diff(ndata->delta.v[2], delta_rot1);

    // We want to treat backward and forward motion symmetrically for the
    // noise model to be applied below.  The standard model seems to assume
    // forward motion.
    delta_rot1_noise = std::min(fabs(angle_diff(delta_rot1,0.0)),
                                fabs(angle_diff(delta_rot1,M_PI)));
    delta_rot2_noise = std::min(fabs(angle_diff(delta_rot2,0.0)),
                                fabs(angle_diff(delta_rot2,M_PI)));

    for (int i = 0; i < set->sample_count; i++)
    {
      pf_sample_t* sample = set->samples + i;

      // Sample pose differences
      delta_rot1_hat = angle_diff(delta_rot1,
                                  pf_ran_gaussian(sqrt(this->alpha1*delta_rot1_noise*delta_rot1_noise +
                                                       this->alpha2*delta_trans*delta_trans)));
      delta_trans_hat = delta_trans - 
              pf_ran_gaussian(sqrt(this->alpha3*delta_trans*delta_trans +
                                   this->alpha4*delta_rot1_noise*delta_rot1_noise +
                                   this->alpha4*delta_rot2_noise*delta_rot2_noise));
      delta_rot2_hat = angle_diff(delta_rot2,
                                  pf_ran_gaussian(sqrt(this->alpha1*delta_rot2_noise*delta_rot2_noise +
                                                       this->alpha2*delta_trans*delta_trans)));

      // Apply sampled update to particle pose
      sample->pose.v[0] += delta_trans_hat * 
              cos(sample->pose.v[2] + delta_rot1_hat);
      sample->pose.v[1] += delta_trans_hat * 
              sin(sample->pose.v[2] + delta_rot1_hat);
      sample->pose.v[2] += delta_rot1_hat + delta_rot2_hat;
    }
  }
  break;
  }
  return true;
}
