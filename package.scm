;;; Copyright (C) 2015 Ricardo Wurmus <ricardo.wurmus@mdc-berlin.de>
;;;
;;; Dorina is free software: you can redistribute it and/or modify it under
;;; the terms of the GNU General Public License as published by the Free
;;; Software Foundation, either version 3 of the License, or (at your option)
;;; any later version.

(use-modules (guix packages)
             ((guix licenses) #:prefix license:)
             (guix build-system python)
             (gnu packages)
             (gnu packages bioinformatics)
             (gnu packages python))

(package
  (name "dorina")
  (version "2.0.0")
  (source #f)
  (build-system python-build-system)
  (arguments
   `(#:python ,python-2
     #:phases
     (modify-phases %standard-phases
       (replace 'check
                (lambda _
                  (zero? (system* "nosetests"
                                  "-v"
                                  "--with-coverage"
                                  "--cover-html"
                                  "--cover-package=\"dorina\"")))))))
  (propagated-inputs
   `(("python2-cython" ,python2-cython)
     ("python2-pybedtools" ,python2-pybedtools)))
  (native-inputs
   `(("python2-coverage" ,python2-coverage)
     ("python2-nose" ,python2-nose)))
  (synopsis "Database of posttranscriptional regulatory elements")
  (description
   "DoRiNA is a database of posttransscriptional regulatory elements.  This
package provides tools to interact with the database.")
  (home-page "https://github.com/dieterich-lab/dorina")
  (license license:gpl3+))
