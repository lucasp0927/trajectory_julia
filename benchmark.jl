#using ProfileView
function benchmark_value(iter::Integer,sfn::ScalarFieldNode{2})
    #TODO add time.
#    Profile.init(delay=0.01)
    info("Benchmarking Fields.value()...")
    srand()
    x_start::Float64 = sfn.position[1]
    y_start::Float64 = sfn.position[2]
    x_size::Float64 = sfn.size[1]
    y_size::Float64 = sfn.size[2]
    t::Float64 = 172
    for i = 1:2
        x = x_start+x_size*rand()
        y = y_start+y_size*rand()
        v = Fields.value([x,y],t,sfn)
    end
#    Profile.clear()
    @time begin
        for i = 1:iter
            x = x_start+x_size*rand()
            y = y_start+y_size*rand()
            v = Fields.value([x,y],t,sfn)
        end
    end
#    ProfileView.view()
#    sleep(100000)
end
