use crate::MAX_DEGREE::MAX_DEGREE;
use crate::MAX_DEGREE::MAX_DEGREE_LOG;
use crate::ROOT_OF_UNITY::ROOT_OF_UNITY;

pub struct Domain {
    elements: [Field; MAX_DEGREE + 1],
}

impl Domain {


    // The root of unity domain is evaluated at compile time
    pub comptime fn new() -> Self {
    let elements = comptime {
        let difference = 28 - MAX_DEGREE_LOG;

        let mut root_of_unity = ROOT_OF_UNITY;
        for _ in 0..difference {
            root_of_unity = root_of_unity * root_of_unity;
        }

        let mut elements = [0; MAX_DEGREE + 1];
        elements[0] = 1;
        elements[1] = root_of_unity;
        for i in 2..MAX_DEGREE + 1 {
            elements[i] = elements[i - 1] * root_of_unity;
        }
        elements
    };
    let result = Self { elements };
    assert(result.check_validity());
    result
    }



    pub fn new_unchecked(elements: [Field; MAX_DEGREE + 1]) -> Self {
        Self { elements }
    }

    pub fn check_validity(self) -> bool {
        assert(self.elements[0] == 1);
        let primitive_root = self.elements[1];
        let mut current = primitive_root;
        for i in 1..MAX_DEGREE + 1 {
            assert(current != 1);
            assert(self.elements[i] == current);

            current = current * primitive_root;

        }

        assert(current == 1);

        true

    }


    pub fn primitive_root(self) -> Field {
        self.elements[1]
    }

    pub fn get_root(self, i: u32) -> Field {
        self.elements[i]
    }

}



//this test only works for MAX_DEGREE = 7
// #[test]
// fn test_domain_size_8() {
//     let domain = Domain::new_unchecked([1,
//     0x2b337de1c8c14f22ec9b9e2f96afef3652627366f8170a0a948dad4ac1bd5e80,
//     0x30644e72e131a029048b6e193fd841045cea24f6fd736bec231204708f703636,
//     0x1d59376149b959ccbd157ac850893a6f07c2d99b3852513ab8d01be8e846a566,
//     0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000000,
//     0x0530d09118705106cbb4a786ead16926d5d174e181a26686af5448492e42a181,
//     0xb3c4d79d41a91758cb49c3517c4604a520cff123608fc9cb,
//     0x130b17119778465cfb3acaee30f81dee20710ead41671f568b11d9ab07b95a9b,
//     ]);


//     assert(domain.check_validity());
// }

#[test]
fn test_new_domain() {
    let domain = comptime { Domain::new() };
    assert(domain.check_validity());
    for i in 0..MAX_DEGREE + 1 {
        println(domain.get_root(i));
    }
}


