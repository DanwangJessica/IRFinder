/* parsestr.c
 *    Routines for parsing and splitting strings at white space.
 *
 * Copyright (C) 2001-2006 The Regents of the University of California
 * by Glenn Tesler
 *
 * See file COPYRIGHT for details.
 *****************************************************************************
 * This file is part of GRIMM-Synteny.
 *
 * GRIMM-Synteny is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License, Version 2,
 * dated June 1991, as published by the Free Software Foundation.
 *
 * GRIMM-Synteny is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

/* Last modified on Sun Sep 3, 2006, by Glenn Tesler
 */

#include <ctype.h>
#include "parsestr.h"

/* count number of space separated fields on line,
 * terminating at EOL or comment ('#')
 */
int count_fields(char *s)
{
  int nf = 0;

  s = skip_blanks(s);
  while (*s != '\0' && *s != '#') {
    nf++;
    s = skip_nonblanks(s);
    s = skip_blanks(s);
  }

  return nf;
}

/* skip over whitespace, return pointer to next char after whitespace */
char *skip_blanks(char *s)
{
  while (*s != '\0' && isspace((int) *s)) s++;
  return s;
}

/* skip over nonwhitespace, return pointer to next char after
 * nonwhitespace or to end of line
 */
char *skip_nonblanks(char *s)
{
  while (*s != '\0' && !isspace((int) *s)) s++;
  return s;
}