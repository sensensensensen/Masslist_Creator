module Masslist_Creator

    # Include the implementation file containing the logic and GUI
    include("Masslist_final.jl")

    # Bring the inner module into scope and export it
    # This allows users to access Masslist_final.main()
    using .Masslist_final
    export Masslist_final

end # module