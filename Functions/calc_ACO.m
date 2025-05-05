function v_ACO = calc_ACO(r, v)

mu = 398600.44189;
h_hat = (cross(r, v) / norm(cross(r, v)));
r_hat = r / norm(r);

v2mag_cir = sqrt(mu/norm(r));

v_ACO_hat = cross(h_hat, r_hat);
v_ACO = v2mag_cir * v_ACO_hat;

end