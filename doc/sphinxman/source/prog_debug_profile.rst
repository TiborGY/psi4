.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2026 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This file is part of Psi4.
.. #
.. # Psi4 is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU Lesser General Public License as published by
.. # the Free Software Foundation, version 3.
.. #
.. # Psi4 is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU Lesser General Public License for more details.
.. #
.. # You should have received a copy of the GNU Lesser General Public License along
.. # with Psi4; if not, write to the Free Software Foundation, Inc.,
.. # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
.. #
.. # @END LICENSE
.. #

.. include:: autodoc_abbr_options_c.rst

.. _`sec:prog_debug_and_profile`:

=======================
Debugging and Profiling
=======================

Debugging
---------

The preferred method for debugging C++ code in Psi4 is with gdb. To run Psi4 in this way, you must specify the Python executable as the program. Arguments are provided using the ``--args`` flag. Here's an example.::
  $~> gdb --args /usr/bin/python3 /path/to/psi4 input.dat

In order to debug properly, Psi4 needs to be built with the ``CMAKE_BUILD_TYPE`` variable set to either ``Debug`` or ``RelWithDebInfo``. These will output debugging symbols that will allow gdb to print line numbers and
inline function names.:

.. code-block:: bash

   > cmake [your options here] -DCMAKE_BUILD_TYPE=Debug
   > cmake [your options here] -DCMAKE_BUILD_TYPE=RelWithDebInfo

Certain symbols may not be output even with this flag set. In general, any template function used should be visible, and anything with the ``PSI_API`` modifier will be visible. Other variables, functions, and classes will
likely be hidden from the user. To make these symbols visible, you must modify a few variables. For an example, see `TiborGY's debug branch <https://github.com/psi4/psi4/compare/master...TiborGY:psi4:toc_dbg>`_.

Also see :ref:`more debugger directions <faq:gdblldb>` 
and a `[presentation] <https://github.com/psi4/PsiCon2020/blob/master/PsiCon2017/Turney-C%2B%2B.pdf>`_ .
If building using ``psi4-path-advisor cmake``, one should run it straight, not within ``eval $(...)``, note the usage command it outputs, then edit the cache file it has produced to change ``CMAKE_BUILD_TYPE`` to ``Debug`` and ``CMAKE_CXX_FLAGS`` to ``-O0``, then execute the noted ``cmake ... -C cache`` command to configure. 

Debug-only code
^^^^^^^^^^^^^^^

To compile C++ that only runs when the build carries debug symbols, branch on the
``psi::have_debug_symbols`` constant from :source:`psi4/include/psi4/pragma.h`. It is ``true`` for
the ``Debug`` and ``RelWithDebInfo`` configurations and ``false`` otherwise.

.. code-block:: cpp

   #include "psi4/pragma.h"

   if constexpr (psi::have_debug_symbols) {
       outfile->Printf("expensive consistency check said %d\n", check());
   }

Prefer ``if constexpr`` over the preprocessor. The discarded branch is still parsed and
type-checked, so the debug code cannot quietly rot as the codebase moves under it, but it is
dropped before code generation and so costs an optimized build nothing. Where the code genuinely
cannot compile in a release configuration, such as an extra data member or an extra include, use
the ``#if PSI4_HAVE_DEBUG_SYMBOLS`` macro instead.

Note this is deliberately not the same test as ``NDEBUG``, which CMake also defines for
``RelWithDebInfo`` even though that configuration does carry symbols. Plain ``assert()``, which
keys off ``NDEBUG``, is therefore inactive in ``RelWithDebInfo``. Out-of-tree consumers such as
plugins are configured separately and build without ``-g`` by default, so the constant is ``false``
there.

VSCode
^^^^^^

When using gdb wath VSCode, you should set the ``program`` entry to the Python executable, just as before. Arguments can then be placed in the ``args`` entry. If you are debugging a C++ plugin or backend code,
the launch type should be ``cppdbg``.
 

Profiling
---------

Instructions on using Psi4 with a profiler, or/and discuss how timers work
in Psi4... they are hierarchical.... I think we have special timers for
parallel blocks.


