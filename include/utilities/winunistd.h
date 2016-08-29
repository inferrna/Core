/****************************************************************************
 *      This library is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public
 *      License as published by the Free Software Foundation; either
 *      version 2.1 of the License, or (at your option) any later version.
 *
 *      This library is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *      Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this library; if not, write to the Free Software
 *      Foundation,Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#ifndef Yafray_win32_unistd_h
#define Yafray_win32_unistd_h


int getopt(int argc, char *argv[], char *opstring) ;
/* static (global) variables that are specified as exported by getopt() */
extern char *optarg;    /* pointer to the start of the option argument  */
extern int   optind;    /* number of the next argv[] to be evaluated    */
extern int   opterr;    /* non-zero if a question mark should be returned  */


#endif
