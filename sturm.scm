#lang racket
; ------------------------------------------------------------------------------
;         Natan DERROITTE - Programation fonctionelle - Projet 1
; ------------------------------------------------------------------------------
(provide p%q)
(provide sturm-chain)
(provide count-roots)
(provide find-roots)

; All polynomes are represented by the list of size (d+1) of their coefficients
; in increasing order of exponent, where d is the highest exponent whose
; coefficient is non-zero.
; For instance, the polynme P(x) = 2 + 3x + x^3  is represented by (2 3 0 1).


; if p and q are respectively the repsentation of P(x) and Q(x),
; (p%q p q) returns the representation of R(x), the remainer of the division
; of P(x) by Q(x).
(define p%q
  (lambda (p q)
    (cond ((is-empty-list? p) '(0))
          ((< (get-degree p) (get-degree q)) (remove-useless-terms p))
          (else
           (let((result (list (/ (highest-term p) (highest-term q)) (- (get-degree p) (get-degree q)))))
             (p%q (sub-poly p (mult-poly-single-term q result)) q))))))


; if p is a polynome, (sturm-chain p) returns the Sturm chain of p by decreasing
; order of powers.
(define sturm-chain
  (lambda (p)
    (chain-aux (list p (deriv p)) p (deriv p))))

; Soit p la liste contenant les i premiers termes de la chaine de Sturm du pylogone P0 tel que p = '(P0, P1, ..., Pi-1),
; a la liste correspondant au i-2ième terme de la chaine de Sturm (ie Pi-2) et b la liste correspondant au i-1ième terme
; de la chaine de Sturm (ie Pi-1), (chain-aux p a b) retourne la chaine de Sturm complète du poylgone P0 par ordre croissant de puissance.
; En particulier, si P' est la dérivé du polygone P, (chain-aux (list P P') P P') donne la chaine de Sturm du polygone P.
(define chain-aux
  (lambda (p a b)
    (let ((remainder (p%q a b)))
          (if (and (zero? (car remainder)) (= 1 (length remainder))) p
              (let ((minus-remainder (map (lambda (x) (* -1 x )) remainder)))
                (chain-aux (append p (list minus-remainder)) b minus-remainder))))))

; if p is a polynome and both a and b are numbers such that a < b,
; (count-roots p a b) returns the number of roots of p
; on ]a b]
(define count-roots
  (lambda (p a b)
    (- (nb-sign-change(evaluate-chain p a)) (nb-sign-change (evaluate-chain p b)))))
 

; if p is a polynome, both a and b are numbers (such that a < b) and eps
; is a positive real, (find-roots p a b eps) returns the ordered list
; of roots of p on the ]a, b] interval with precision eps
(define find-roots
  (lambda (p a b eps)
    (let ((nbRoots (count-roots p a b))(c (/ (+ a b)2)))
      (cond ((zero? nbRoots) '())
            ((= 1 nbRoots) (list (find-root p a b eps)))
            ((< (- b a) eps) (list c))
            ((zero? (evaluate-poly p c)) (cons (find-roots p a (+ c eps) eps) (find-roots p (+ c eps) b eps)))
            (else (append (find-roots p a c eps) (find-roots p c b eps)))))))
        

; Soit a et b deux réels, p un polygone P(x) n'ayant qu'une racine dans l'interval ]a , b[ et eps une précision, (find-root p a b eps) a pour valeur la racine de
; P(x) à une précision eps près. Cette fonction utilise une méthode dichotomique.
(define find-root
  (lambda (p a b eps)
    (let ((c (/ (+ a b) 2))(fc (evaluate-poly p (/ (+ a b) 2))))
      (cond ((or(zero? fc)(< (/ (- b a) 2) eps)) c)
            ((= (sgn fc) (sgn (evaluate-poly p a))) (find-root p c b eps))
            (else (find-root p a c eps))))))

; Soit l une liste quelconque, le prédicat (is-empty-list? l) a pour valeur vrai si l est la liste nulle ou uniquement composée de 0 et faux sinon.
(define is-empty-list?
  (lambda (l)
    (cond
      ((null? l) #t)
      ((not (zero? (car l ))) #f)
      (else (is-empty-list? (cdr l))))))
        
; Soit l une liste quelconque (get-degree l) a pour valeur l'index du dernier élément non nul de la liste l, commençant a compter le premier élément de l à l'index 0.
; En particulier, si l représente une poylgone, (get-degree l) a pour valeur le
; degré du polygone.
(define get-degree
  (lambda (l)
    (if (is-empty-list? l) 0
        (get-degree-aux l -1))))

; Soit l une liste quelconque et n un entier, (get-degree-aux l n) a pour valeur n + l'index du dernier élément non nul de l,
; commençant a compter le premier élément de l  à l'index 0.
; En particulier, si l ne contient aucun non-nul  (get-degree-aux l n) à pour valeur n.
(define get-degree-aux
  (lambda (l n)
    (cond
      ((is-empty-list? l) n)
      (else (get-degree-aux (cdr l) (+ n 1))))))

; Soit l une liste quelconque (highest-term l) a pour valeur le dernier élément non nul de la liste l, commençant a compter le premier élément de l à l'index 0.
; En particulier, si l représente une poylgone, (highest-term l) a pour valeur le coefficient du terme de plus haut degré du polygone.
(define highest-term
  (lambda (l)
    (cond
      ((null? l) '()) ;Should never happen.
      ((is-empty-list? (cdr l)) (car l))
      (else (highest-term (cdr l))))))

; Soit p et q,deux polynomes P(x) et Q(x) quelconques, (sub-poly p q) retourne la liste correspndant au polygone resultant de la soustraction de P et Q
(define sub-poly
  (lambda (p q)
    (cond ((= (length p)(length q)) (map - p q))
          ((< (length p)(length q)) (sub-poly (add-zeros-terms-end p (- (length q) (length p))) q))
          (else (sub-poly p (add-zeros-terms-end q (- (length p) (length q))))))))

; Soit la liste p quelconque et un entier n quelconque, (add-zeros-terms-end p n) retourne la liste p concaténée de n 0
(define add-zeros-terms-end
  (lambda (p n)
    (if (zero? n) p
        (add-zeros-terms-end (append p '(0)) (- n 1)))))

; Soit la liste p quelconque et un entier n quelconque, (add-zeros-terms-end p n) retourne la liste comprenant n 0 suivie de p
(define add-zeros-terms-begining
  (lambda (p n)
    (if (zero? n) p
        (add-zeros-terms-begining (append '(0) p) (- n 1)))))

; Soit p un polygone P(x) quelconque  et q une liste de deux éléments correspondant  à un terme en x, q(x), où le premier élément est son coefficient et le deuxième son degrée,
; (mult-poly-single-term p q) retourne la liste correspondante au polygone résultant de la multiplication de P(x) par le terme q(x)
(define mult-poly-single-term
  (lambda (p q)
    (map (lambda (x) (* x (car q))) (add-zeros-terms-begining p (cadr q)))))

; Soit p une liste quelconque, (remove-useless-terms p) a pour valeur la liste p où tous les  éventuels 0 terminant cette liste ont été retirés
(define remove-useless-terms
  (lambda (p)
    (if (is-empty-list? p) '()
        (cons (car p) (remove-useless-terms (cdr p))))))

; Soit p un polygone P(x) quelconque (deriv p) a pour valeur la liste au polygone P'(x) où  P'(x) est la dérviée de P(x) selon x.
(define deriv
  (lambda (p)
    (if (null? p) p
        (deriv-aux 1 (cdr p)))))

; Soit n un entier et p un polygone P(x) quelconque dont le premier élément est de degré d (et non 0), (deriv-aux n p) correspondant à la dérivée de P(x), ie P'(x)
(define deriv-aux
  (lambda (n p)
    (if (null? p) '()
        (cons (* n (car p)) (deriv-aux (+ n 1) (cdr p))))))

; Soit p un polygone P(x) quelconque et a un réel quelconque, (evaluate-poly p a) a pour valeur P(a)
(define evaluate-poly
  (lambda (p a)
    (evaluate-poly-aux p a 0 0)))

; Soit p un polygone P(x) quelconque dont le premier terme est de degrée d (et non 0), a et r deux réels quelconques,  d un entier correspondant au degrée en question,
; (evaluate-poly-aux p a r d) a pour valeur r + P(a)
(define evaluate-poly-aux
  (lambda ( p a r d)
    (if (null? p) r
        (evaluate-poly-aux (cdr p) a (+ r (* (car p) (expt a d))) (+ 1 d)))))

; Soit p un polygone P(x) quelconque et a un réel quelconque, (evaluate-chain p a) a pour valeur à la liste correspondant à la série de Sturm de P(x) évaluée en a
(define evaluate-chain
  (lambda (p a)
    (evaluate-chain-aux (sturm-chain p) a '())))

; Soit chain une liste de poylgone, a un réel quelconque et l une liste quelconque, (evaluate-chain-aux chain a l) a pour valeur la liste l concaténée des évaluations des
; polygones contenus dans chain en a. Les évaluations sont concaténées dans l'ordre où les polygones sont contenus dans chain
(define evaluate-chain-aux
  (lambda (chain a l)
    (if (null? chain) l
        (evaluate-chain-aux (cdr chain) a (append l (list(evaluate-poly (car chain) a)))))))

;Soit l une liste quelconque contenant des nombres, (nb-sign-change l) a pour valeur le nombre de changement de signes des éléments de l.
(define nb-sign-change
  (lambda (l)
    (nb-sign-change-aux (cdr l) 0 (sgn (car l)))))

; Soit l une liste quelconque contenant des nombres, (nb-sign-change l) a pour valeur n + (le nombre de changement de signes des éléments de l) + 1 si le premier éléments de l
; n'est pas de même signe que s.
(define nb-sign-change-aux
  (lambda (l n s)
    (cond ((null? l) n)
          ((or(zero? (car l ))(= (sgn(car l)) s)) (nb-sign-change-aux (cdr l) n s))
          (else (nb-sign-change-aux (cdr l) (+ n 1) (sgn (car l)))))))