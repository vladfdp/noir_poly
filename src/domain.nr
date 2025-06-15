use crate::ROOT_OF_UNITY::ROOT_OF_UNITY;

pub struct Domain<let DOMAIN_SIZE: u32, let DOMAIN_LOG_SIZE: u32> {
    elements: [Field; DOMAIN_SIZE],
}

impl<let DOMAIN_SIZE: u32, let DOMAIN_LOG_SIZE: u32> Eq for Domain<DOMAIN_SIZE, DOMAIN_LOG_SIZE> {
    fn eq(self, other: Self) -> bool {
        (self.elements == other.elements)
    }
}

impl<let DOMAIN_SIZE: u32, let DOMAIN_LOG_SIZE: u32> Domain<DOMAIN_SIZE, DOMAIN_LOG_SIZE> {


    // The root of unity domain is evaluated at compile time
    pub comptime fn new() -> Self {
        assert(DOMAIN_LOG_SIZE <= 28);
        assert(DOMAIN_SIZE == 1 << (DOMAIN_LOG_SIZE as u8));

        let difference = 28 - DOMAIN_LOG_SIZE;

        let mut root_of_unity = ROOT_OF_UNITY;
        for _ in 0..difference {
            root_of_unity = root_of_unity * root_of_unity;
        }

        let mut elements = [0; DOMAIN_SIZE];
        elements[0] = 1;
        elements[1] = root_of_unity;
        for i in 2..DOMAIN_SIZE {
            elements[i] = elements[i - 1] * root_of_unity;
        }
        let result = Self { elements };
        result
    }

 



    pub fn new_unchecked(elements: [Field; DOMAIN_SIZE]) -> Self {
        Self { elements }
    }

    pub fn check_validity(self) -> bool {
        assert(self.elements[0] == 1);
        let primitive_root = self.elements[1];
        let mut current = primitive_root;
        for i in 1..DOMAIN_SIZE {
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


#[test]
fn test_domain_size_8() {
    let domain = Domain::<8, 3>::new_unchecked([1,
    0x2b337de1c8c14f22ec9b9e2f96afef3652627366f8170a0a948dad4ac1bd5e80,
    0x30644e72e131a029048b6e193fd841045cea24f6fd736bec231204708f703636,
    0x1d59376149b959ccbd157ac850893a6f07c2d99b3852513ab8d01be8e846a566,
    0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000000,
    0x0530d09118705106cbb4a786ead16926d5d174e181a26686af5448492e42a181,
    0xb3c4d79d41a91758cb49c3517c4604a520cff123608fc9cb,
    0x130b17119778465cfb3acaee30f81dee20710ead41671f568b11d9ab07b95a9b,
    ]);

    let domain2 = comptime { Domain::<8, 3>::new() };

    assert(domain == domain2);
}

#[test]
fn test_new_domain() {
    let domain = comptime { Domain::<256, 8>::new() };
    assert(domain.check_validity());
}

#[test]
fn test_primitive_root() {
    let domain = comptime { Domain::<8, 3>::new() };
    let primitive_root = domain.primitive_root();
    
    // Verify primitive root property: root^8 = 1 and root^k != 1 for k < 8
    let mut power = primitive_root;
    for i in 1..8 {
        assert(power != 1);
        power = power * primitive_root;
    }
    assert(power == 1);
}

