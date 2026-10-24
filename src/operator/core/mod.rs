// helpers; not in prelude, called by path
pub mod flow_helper;

// scalar domain operators
pub mod opscl_dom_adv;
pub mod opscl_dom_advunity;
pub mod opscl_dom_diff;
pub mod opscl_dom_div;
pub mod opscl_dom_pref;
pub mod opscl_dom_prod;
pub mod opscl_dom_pstab;
pub mod opscl_dom_src;
pub mod opscl_dom_time;
pub mod opscl_dom_timeprev;
pub mod opscl_dom_timeunity;

// scalar boundary operators
pub mod opscl_bnd_adv;
pub mod opscl_bnd_advunity;
pub mod opscl_bnd_diff;
pub mod opscl_bnd_dir;
pub mod opscl_bnd_div;
pub mod opscl_bnd_neu;
pub mod opscl_bnd_prod;
pub mod opscl_bnd_pstab;
pub mod opscl_bnd_src;
pub mod opscl_bnd_time;
pub mod opscl_bnd_timeprev;
pub mod opscl_bnd_timeunity;

// scalar interface operators
pub mod opscl_itf_adv;
pub mod opscl_itf_advunity;
pub mod opscl_itf_cont;
pub mod opscl_itf_diff;
pub mod opscl_itf_div;
pub mod opscl_itf_prod;
pub mod opscl_itf_pstab;
pub mod opscl_itf_src;
pub mod opscl_itf_time;
pub mod opscl_itf_timeprev;
pub mod opscl_itf_timeunity;

// vector domain operators
pub mod opvec_dom_adv;
pub mod opvec_dom_ddil;
pub mod opvec_dom_diff;
pub mod opvec_dom_dtrn;
pub mod opvec_dom_grad;
pub mod opvec_dom_src;
pub mod opvec_dom_time;

// vector boundary operators
pub mod opvec_bnd_adv;
pub mod opvec_bnd_ddil;
pub mod opvec_bnd_diff;
pub mod opvec_bnd_dir;
pub mod opvec_bnd_dtrn;
pub mod opvec_bnd_grad;
pub mod opvec_bnd_trac;
pub mod opvec_bnd_src;
pub mod opvec_bnd_time;

// vector interface operators
pub mod opvec_itf_adv;
pub mod opvec_itf_cont;
pub mod opvec_itf_ddil;
pub mod opvec_itf_diff;
pub mod opvec_itf_dtrn;
pub mod opvec_itf_grad;
pub mod opvec_itf_src;
pub mod opvec_itf_time;
