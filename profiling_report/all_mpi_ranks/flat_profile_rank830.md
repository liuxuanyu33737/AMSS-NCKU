# MPI 进程 830 的 gprof 性能报告

- 数据文件：`gmon.out.830`
- 函数数量：129
- 分析工具：GNU gprof

|排名|占比(%)|累计占比(%)|自身时间(s)|调用次数|函数名称|
|---:|---:|---:|---:|---:|---|
|1|57.50|57.50|223.33|5457|`compute_rhs_bssn_`|
|2|8.98|66.48|34.86|130968|`lopsided_`|
|3|7.15|73.63|27.78|60027|`fdderivs_`|
|4|5.95|79.58|23.11|130968|`kodis_`|
|5|5.40|84.98|20.97|91237|`fderivs_`|
|6|5.05|90.03|19.63|207686|`prolong3_`|
|7|2.25|92.28|8.72|5280|`enforce_ga_`|
|8|2.20|94.48|8.53|126720|`rungekutta4_rout_`|
|9|1.75|96.23|6.81|752876|`symmetry_bd_`|
|10|0.85|97.08|3.30|130070|`restrict3_`|
|11|0.82|97.90|3.19||`_init`|
|12|0.50|98.40|1.96|15936|`average2_`|
|13|0.32|98.72|1.25|521994240|`misc::fact(int)`|
|14|0.23|98.95|0.91|240|`admmass_bssn_`|
|15|0.20|99.15|0.79|1579063|`copy_`|
|16|0.20|99.35|0.78|20|`getnp4_`|
|17|0.08|99.43|0.31|46448640|`misc::Wigner_d_function(int, int, int, double)`|
|18|0.08|99.51|0.30|1760|`Patch::Interp_Points(MyList<var>*, int, double**, double*, int)`|
|19|0.07|99.58|0.28|566136|`Ansorg::interpolate_tri_bar(double, double, double, int, int, int, double*, double*, double*, double*)`|
|20|0.06|99.64|0.24|1802560|`polint_`|
|21|0.06|99.70|0.24|5280|`lowerboundset_`|
|22|0.06|99.76|0.23|240|`surface_integral::surf_Wave(double, int, cgh*, var*, var*, int, int, int, double*, double*, monitor*)`|
|23|0.03|99.79|0.12|870502|`monitor::monitor(char const*, int, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)`|
|24|0.03|99.82|0.12|1320|`bssn_class::Step(int, int)`|
|25|0.03|99.85|0.11|124800|`sommerfeld_rout_`|
|26|0.03|99.88|0.11|29043|`Parallel::transfer(MyList<Parallel::gridseg>**, MyList<Parallel::gridseg>**, MyList<var>*, MyList<var>*, int)`|
|27|0.02|99.90|0.09|1920|`sommerfeld_routbam_`|
|28|0.02|99.92|0.06|240|`surface_integral::surf_MassPAng(double, int, cgh*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, double*, monitor*)`|
|29|0.01|99.93|0.05|1260|`l2normhelper_`|
|30|0.01|99.94|0.03|13214167|`Block::getdX(int)`|
|31|0.01|99.95|0.03|9|`get_ansorg_nbhs_`|
|32|0.01|99.96|0.02|230517|`Patch::Interp_ONE_Point(MyList<var>*, double*, double*, int)`|
|33|0.01|99.97|0.02|150168|`Parallel::gs_subtract(MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>*)`|
|34|0.01|99.98|0.02|120|`average_`|
|35|0.00|99.98|0.01|566136|`Ansorg::ps_u_at_xyz(double, double, double)`|
|36|0.00|99.98|0.01|140352|`Parallel::build_owned_gsl2(Patch*, int)`|
|37|0.00|99.98|0.01|41920|`decide3d_`|
|38|0.00|99.98|0.01|35552|`Parallel::build_owned_gsl4(Patch*, int, int)`|
|39|0.00|99.98|0.01|14923|`Parallel::Sync(Patch*, MyList<var>*, int)`|
|40|0.00|99.98|0.01|8212|`Parallel::Sync(MyList<Patch>*, MyList<var>*, int)`|
|41|0.00|99.98|0.01|1093|`Parallel::prepare_inter_time_level(Patch*, MyList<var>*, MyList<var>*, MyList<var>*, MyList<var>*, int)`|
|42|0.00|99.98|0.00|1698408|`Ansorg::find_point_bisection(double, int, double*, int)`|
|43|0.00|99.98|0.00|701658|`Patch::getdX(int)`|
|44|0.00|99.98|0.00|566136|`Ansorg::xyz_to_ABp(double, double, double, double*, double*, double*)`|
|45|0.00|99.98|0.00|231952|`Parallel::build_gstl(MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>**, MyList<Parallel::gridseg>**)`|
|46|0.00|99.98|0.00|120144|`Parallel::build_owned_gsl0(Patch*, int)`|
|47|0.00|99.98|0.00|65696|`Parallel::build_bulk_gsl(Block*, Patch*)`|
|48|0.00|99.98|0.00|42240|`Block::swapList(MyList<var>*, MyList<var>*, int)`|
|49|0.00|99.98|0.00|41920|`global_interp_`|
|50|0.00|99.98|0.00|41920|`polin3_`|
|51|0.00|99.98|0.00|21821|`Parallel::build_complete_gsl(Patch*)`|
|52|0.00|99.98|0.00|21000|`cgh::Interp_One_Point(MyList<var>*, double*, double*, int)`|
|53|0.00|99.98|0.00|19327|`Parallel::gsl_subtract(MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>*)`|
|54|0.00|99.98|0.00|19278|`Parallel::build_buffer_gsl(Patch*)`|
|55|0.00|99.98|0.00|14923|`Parallel::build_ghost_gsl(Patch*)`|
|56|0.00|99.98|0.00|4355|`Parallel::OutBdLow2Hi(Patch*, Patch*, MyList<var>*, MyList<var>*, int)`|
|57|0.00|99.98|0.00|1939|`misc::parse_parts(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, int&)`|
|58|0.00|99.98|0.00|1406|`Parallel::Restrict(MyList<Patch>*, MyList<Patch>*, MyList<var>*, MyList<var>*, int)`|
|59|0.00|99.98|0.00|1321|`setpbh(int, double**, double*, int)`|
|60|0.00|99.98|0.00|1308|`bssn_class::Constraint_Out()`|
|61|0.00|99.98|0.00|1280|`Parallel::PatList_Interp_Points(MyList<Patch>*, MyList<var>*, int, double**, double*, int)`|
|62|0.00|99.98|0.00|1260|`Parallel::L2Norm(Patch*, var*)`|
|63|0.00|99.98|0.00|700|`cgh::Regrid_Onelevel(int, int, int, double**, double**, MyList<var>*, MyList<var>*, MyList<var>*, MyList<var>*, bool, monitor*)`|
|64|0.00|99.98|0.00|640|`bssn_class::compute_Porg_rhs(double**, double**, var*, var*, var*, int)`|
|65|0.00|99.98|0.00|639|`misc::parse_parts(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, int&, int&, int&)`|
|66|0.00|99.98|0.00|464|`Block::Block(int, int*, double*, int, int, int, int, int)`|
|67|0.00|99.98|0.00|464|`Block::~Block()`|
|68|0.00|99.98|0.00|440|`monitor::writefile(double, int, double*)`|
|69|0.00|99.98|0.00|240|`monitor::writefile(double, int, double*, double*)`|
|70|0.00|99.98|0.00|167|`var::var(char const*, int, double, double, double)`|
|71|0.00|99.98|0.00|167|`var::~var()`|
|72|0.00|99.98|0.00|107|`Patch::Patch(int, int*, double*, int, bool, int)`|
|73|0.00|99.98|0.00|107|`Patch::~Patch()`|
|74|0.00|99.98|0.00|107|`Parallel::partition3(int*, int, int*, int, int*)`|
|75|0.00|99.98|0.00|58|`cgh::construct_patchlist(int, int)`|
|76|0.00|99.98|0.00|58|`Parallel::KillBlocks(MyList<Patch>*)`|
|77|0.00|99.98|0.00|58|`Parallel::distribute(MyList<Patch>*, int, int, int, bool, int)`|
|78|0.00|99.98|0.00|58|`Parallel::add_ghost_touch(MyList<Parallel::gridseg>*&)`|
|79|0.00|99.98|0.00|58|`Parallel::cut_gsl(MyList<Parallel::gridseg>*&)`|
|80|0.00|99.98|0.00|58|`Parallel::merge_gsl(MyList<Parallel::gridseg>*&, double)`|
|81|0.00|99.98|0.00|53|`Parallel::checkgsl(MyList<Parallel::gridseg>*, bool)`|
|82|0.00|99.98|0.00|49|`cgh::recompose_cgh_Onelevel(int, int, MyList<var>*, MyList<var>*, MyList<var>*, MyList<var>*, int, bool)`|
|83|0.00|99.98|0.00|49|`Parallel::fill_level_data(MyList<Patch>*, MyList<Patch>*, MyList<Patch>*, MyList<var>*, MyList<var>*, MyList<var>*, MyList<var>*, int, bool, bool)`|
|84|0.00|99.98|0.00|49|`Parallel::build_complete_gsl_virtual(MyList<Patch>*)`|
|85|0.00|99.98|0.00|49|`Parallel::cut_gs(MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>*&)`|
|86|0.00|99.98|0.00|21|`bssn_class::Interp_Constraint(bool)`|
|87|0.00|99.98|0.00|21|`bssn_class::Compute_Constraint()`|
|88|0.00|99.98|0.00|21|`cgh::cgh(int, int, int, char*, int, monitor*)`|
|89|0.00|99.98|0.00|20|`bssn_class::Compute_Psi4(int)`|
|90|0.00|99.98|0.00|20|`bssn_class::RecursiveStep(int)`|
|91|0.00|99.98|0.00|20|`perf::MemoryUsage(unsigned long*, unsigned long*, unsigned long*, unsigned long*, unsigned long*, unsigned long*, int)`|
|92|0.00|99.98|0.00|20|`perf::sample_mem_usage(int)`|
|93|0.00|99.98|0.00|15|`var::setpropspeed(double)`|
|94|0.00|99.98|0.00|13|`Parallel::aligncheck(double*, double*, int, double*, int*)`|
|95|0.00|99.98|0.00|12|`Patch::checkPatch(bool)`|
|96|0.00|99.98|0.00|12|`Parallel::Dump_Data(Patch*, MyList<var>*, char*, double, double, int)`|
|97|0.00|99.98|0.00|9|`Parallel::Dump_Data(MyList<Patch>*, MyList<var>*, char*, double, double)`|
|98|0.00|99.98|0.00|9|`std::_Rb_tree<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::_Select1st<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > >::find(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&)`|
|99|0.00|99.98|0.00|7|`Parallel::prepare_inter_time_level(Patch*, MyList<var>*, MyList<var>*, MyList<var>*, int)`|
|100|0.00|99.98|0.00|5|`monitor::~monitor()`|
|101|0.00|99.98|0.00|3|`frame_dummy`|
|102|0.00|99.98|0.00|2|`checkpoint::addvariablelist(MyList<var>*)`|
|103|0.00|99.98|0.00|2|`misc::inversearray(double*, int)`|
|104|0.00|99.98|0.00|2|`std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >::~pair()`|
|105|0.00|99.98|0.00|2|`std::pair<std::_Rb_tree_iterator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, bool> std::_Rb_tree<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::_Select1st<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > >::_M_insert_unique<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >(std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >&&)`|
|106|0.00|99.98|0.00|1|`bssn_class::Initialize()`|
|107|0.00|99.98|0.00|1|`bssn_class::Read_Ansorg()`|
|108|0.00|99.98|0.00|1|`bssn_class::Setup_Black_Hole_position()`|
|109|0.00|99.98|0.00|1|`bssn_class::Evolve(int)`|
|110|0.00|99.98|0.00|1|`bssn_class::bssn_class(double, double, double, double, double, double, double, int, int, char*, double, double, double, int, int, int, double, double)`|
|111|0.00|99.98|0.00|1|`bssn_class::~bssn_class()`|
|112|0.00|99.98|0.00|1|`checkpoint::checkpoint(bool, char const*, int)`|
|113|0.00|99.98|0.00|1|`checkpoint::~checkpoint()`|
|114|0.00|99.98|0.00|1|`surface_integral::surface_integral(int)`|
|115|0.00|99.98|0.00|1|`surface_integral::~surface_integral()`|
|116|0.00|99.98|0.00|1|`cgh::compose_cgh(int)`|
|117|0.00|99.98|0.00|1|`cgh::settrfls(int)`|
|118|0.00|99.98|0.00|1|`cgh::read_bbox(int, char*)`|
|119|0.00|99.98|0.00|1|`cgh::sethandle(monitor*)`|
|120|0.00|99.98|0.00|1|`cgh::~cgh()`|
|121|0.00|99.98|0.00|1|`misc::gaulegf(double, double, double*, double*, int)`|
|122|0.00|99.98|0.00|1|`perf::perf()`|
|123|0.00|99.98|0.00|1|`perf::~perf()`|
|124|0.00|99.98|0.00|1|`Ansorg::set_ABp()`|
|125|0.00|99.98|0.00|1|`Ansorg::Ansorg(char*, int)`|
|126|0.00|99.98|0.00|1|`Ansorg::~Ansorg()`|
|127|0.00|99.98|0.00|1|`monitor::print_message(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)`|
|128|0.00|99.98|0.00|1|`std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, int, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, int> > >::~map()`|
|129|0.00|99.98|0.00|1|`std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >::pair<char const (&) [9], std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, true>(char const (&) [9], std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&)`|
